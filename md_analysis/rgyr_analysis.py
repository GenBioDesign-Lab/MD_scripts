import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Get input files from command line
if len(sys.argv) == 3:
    psf_file = sys.argv[1]
    traj_file = sys.argv[2]
else:
    print("Usage: python rgyr_analysis.py <path_to_psf_file> <path_to_trajectory_file>")
    sys.exit(1)

# Check if files exist
if not os.path.exists(psf_file):
    print(f"Error: PSF file {psf_file} not found")
    sys.exit(1)
if not os.path.exists(traj_file):
    print(f"Error: Trajectory file {traj_file} not found")
    sys.exit(1)

# Create output directory
traj_dir = os.path.dirname(traj_file)
output_dir = os.path.join(traj_dir, "output")
os.makedirs(output_dir, exist_ok=True)

# Get base name of trajectory file for output naming
traj_base = os.path.splitext(os.path.basename(traj_file))[0]

# Load data
u = mda.Universe(psf_file, traj_file)
protein = u.select_atoms("protein") #if you want choose another atom, just change the name, prefer to this link: https://userguide.mdanalysis.org/stable/selections.html

rgyr = []
times = []

for ts in u.trajectory:
    rgyr.append(protein.radius_of_gyration())
    times.append(ts.time)

df = pd.DataFrame({
    "Time (ps)": times,
    "Radius of Gyration (Å)": rgyr
})

# Save data to CSV in output directory
csv_path = os.path.join(output_dir, f"{traj_base}_rgyr.csv")
df.to_csv(csv_path, index=False)

# Create plot
plt.figure(figsize=(10, 6))
plt.plot(times, rgyr)
plt.xlabel("Time (ps)")
plt.ylabel("Radius of Gyration (Å)")
plt.title(f"Protein Compactness - {traj_base}")

# Save plot to output directory
plot_path = os.path.join(output_dir, f"{traj_base}_rgyr_plot.png")
plt.savefig(plot_path, dpi=300)

# Show the plot
plt.show()