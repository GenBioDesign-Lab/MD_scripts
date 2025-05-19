import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Get input files either from command line or user input
if len(sys.argv) == 3:
    psf_file = sys.argv[1]
    dcd_file = sys.argv[2]
else:
    print("Usage: python rmsd_analysis_manual.py <path_to_psf_file> <path_to_dcd_file>")
    print("Please enter the paths manually:")
    psf_file = input("Enter path to PSF file: ")
    dcd_file = input("Enter path to DCD file: ")

# Check if files exist
if not os.path.exists(psf_file):
    print(f"Error: PSF file {psf_file} not found")
    sys.exit(1)
if not os.path.exists(dcd_file):
    print(f"Error: DCD file {dcd_file} not found")
    sys.exit(1)

# Load trajectory
u = mda.Universe(psf_file, dcd_file)

# Select atoms for alignment and RMSD calculation (CA atoms)
ref = u.select_atoms("name CA")
mobile = u.select_atoms("name CA")

# Store first frame as reference
u.trajectory[0]
ref_positions = ref.positions.copy()

# Calculate RMSD for each frame
times = []
rmsd_values = []

for ts in u.trajectory:
    # Get current frame time
    times.append(ts.time)
    
    # Calculate RMSD after optimal alignment
    # This internally does the alignment and RMSD calculation in one step
    rmsd_value = mda.analysis.rms.rmsd(mobile.positions, ref_positions, 
                                      center=True, superposition=True)
    rmsd_values.append(rmsd_value)

# Create dataframe and save results
df = pd.DataFrame({
    "Time (ps)": times,
    "RMSD (Å)": rmsd_values
})
df.to_csv("rmsd_results_manual.csv", index=False)

# Plot
plt.plot(df["Time (ps)"], df["RMSD (Å)"])
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (Å)")
plt.title("Backbone RMSD (Manual Calculation)")
plt.savefig("rmsd_plot_manual.png")

# Show the plot
plt.show() 