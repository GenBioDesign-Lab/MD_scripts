import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import pca
import sys
import os

#Usage: python pca_plots.py <psf_file> <trajectory_file>
#Original script: https://userguide.mdanalysis.org/stable/examples/analysis/reduced_dimensions/pca.html

# Check command line arguments
if len(sys.argv) != 3:
    print("Usage: python pca_plots.py <psf_file> <trajectory_file>")
    sys.exit(1)

psf_file = sys.argv[1]
trajectory_file = sys.argv[2]

# Check if files exist
if not os.path.exists(psf_file):
    print(f"Error: PSF file '{psf_file}' not found")
    sys.exit(1)
if not os.path.exists(trajectory_file):
    print(f"Error: Trajectory file '{trajectory_file}' not found")
    sys.exit(1)

# Create output directory in the same folder as the trajectory file
trajectory_dir = os.path.dirname(os.path.abspath(trajectory_file))
output_dir = os.path.join(trajectory_dir, 'output')
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory created at: {output_dir}")

print(f"Running PCA analysis on:\nPSF: {psf_file}\nTrajectory: {trajectory_file}")

# Load trajectory
u = mda.Universe(psf_file, trajectory_file)

# Select backbone atoms
backbone = u.select_atoms('resname CNT')
n_bb = len(backbone)
print(f'There are {n_bb} backbone atoms in the analysis')

# Run PCA on backbone atoms
pc = pca.PCA(u, select='resname CNT', align=True, mean=None, n_components=None)
pc.run()

# Print PCA components shape
print(pc.results.p_components.shape)

# Print first PC variance
print(f"PC1: {pc.results.variance[0]:.5f}")

# Print cumulative variance for first 3 components
for i in range(3):
    print(f"Cumulated variance: {pc.results.cumulated_variance[i]:.3f}")

# Get base name of trajectory file for output naming
base_name = os.path.basename(trajectory_file).split('.')[0]

# Plot cumulative variance
plt.figure(figsize=(10, 6))
plt.plot(pc.results.cumulated_variance[:10])
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance')
plt.title('Cumulative Variance by Principal Component')
plt.savefig(os.path.join(output_dir, f'{base_name}_cumulative_variance.png'), dpi=300, bbox_inches='tight')
plt.close()

# Transform to get projections onto PCs - now using 5 components if available
n_components = min(5, len(pc.results.p_components))
print(f"Using {n_components} components for transformation")
transformed = pc.transform(backbone, n_components=n_components)
print(f"Transformed data shape: {transformed.shape}")

# Create DataFrame for easier analysis
# Only use first 3 for visualization regardless of how many we computed
viz_components = min(3, n_components)
df = pd.DataFrame(transformed[:, :viz_components], 
                 columns=[f'PC{i+1}' for i in range(viz_components)])
df['Time (ps)'] = df.index * u.trajectory.dt
print(df.head())

# Plot PC scatter plots using seaborn PairGrid
g = sns.PairGrid(df, hue='Time (ps)', palette=sns.color_palette('viridis', n_colors=len(df)))
g.map(plt.scatter, marker='.')
plt.savefig(os.path.join(output_dir, f'{base_name}_pca_pairgrid.png'), dpi=300, bbox_inches='tight')
plt.close()

# Plot first two principal components
plt.figure(figsize=(10, 8))
plt.scatter(df['PC1'], df['PC2'], c=df['Time (ps)'], cmap='viridis')
plt.colorbar(label='Time (ps)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title(f'PCA of {base_name} Backbone')
plt.savefig(os.path.join(output_dir, f'{base_name}_pc1_pc2.png'), dpi=300, bbox_inches='tight')
plt.close()

# Plot PC1 vs time
plt.figure(figsize=(12, 6))
plt.plot(df['Time (ps)'], df['PC1'])
plt.xlabel('Time (ps)')
plt.ylabel('PC1')
plt.title('PC1 vs Time')
plt.savefig(os.path.join(output_dir, f'{base_name}_pc1_vs_time.png'), dpi=300, bbox_inches='tight')
plt.close()

# Project the trajectory onto the first principal component
pc1 = pc.results.p_components[:, 0]
trans1 = transformed[:, 0]
projected = np.outer(trans1, pc1) + pc.mean.flatten()
coordinates = projected.reshape(len(trans1), -1, 3)

# Save dataframe to CSV for further analysis
df.to_csv(os.path.join(output_dir, f'{base_name}_pca_data.csv'), index=False)

# Add cosine content analysis
print("\n=== Measuring convergence with cosine content ===")
print("Cosine content measures similarity to a cosine shape (0-1)")
print("Values close to 1 may indicate poor sampling or random diffusion")
print("Values below 0.7 do not indicate poor sampling")

# Calculate cosine content for available PCs (use actual number of components)
available_components = transformed.shape[1]
print(f"Calculating cosine content for {available_components} components")
cosine_contents = []
for i in range(available_components):
    cc = pca.cosine_content(transformed, i)
    cosine_contents.append(cc)
    print(f"Cosine content for PC {i+1} = {cc:.3f}")

# Create a dataframe with cosine content values
cc_df = pd.DataFrame({
    'PC': [f'PC{i+1}' for i in range(available_components)],
    'Cosine Content': cosine_contents
})
cc_df.to_csv(os.path.join(output_dir, f'{base_name}_cosine_content.csv'), index=False)

# Plot cosine content as a bar chart
plt.figure(figsize=(10, 6))
plt.bar(cc_df['PC'], cc_df['Cosine Content'])
plt.axhline(y=0.7, color='r', linestyle='--', label='0.7 threshold')
plt.xlabel('Principal Component')
plt.ylabel('Cosine Content')
plt.title('Cosine Content of Principal Components')
plt.legend()
plt.savefig(os.path.join(output_dir, f'{base_name}_cosine_content.png'), dpi=300, bbox_inches='tight')
plt.close()

# Plot transformed components over time
# Melt the dataframe into a tidy format for plotting
melted = pd.melt(df, id_vars=['Time (ps)'], 
                 value_vars=[f'PC{i+1}' for i in range(viz_components)],
                 var_name="PC", value_name="Value")

# Create FacetGrid to plot PC trajectories
plt.figure(figsize=(12, 8))
g = sns.FacetGrid(melted, col='PC', height=4, aspect=1.2)
g.map_dataframe(sns.lineplot, x='Time (ps)', y='Value')
g.set_axis_labels('Time (ps)', 'Value')
g.fig.suptitle('PC Trajectories Over Time', y=1.05)
plt.savefig(os.path.join(output_dir, f'{base_name}_pc_trajectories.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"\nPCA analysis completed. All results saved to {output_dir}/")
