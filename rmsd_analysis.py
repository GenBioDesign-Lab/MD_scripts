import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Get input files either from command line or user input
if len(sys.argv) == 3:
    psf_file = sys.argv[1]
    dcd_file = sys.argv[2]
else:
    print("Usage: python rmsd_analysis.py <path_to_psf_file> <path_to_dcd_file>")
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

u = mda.Universe(psf_file, dcd_file)

rmsd = RMSD(u, select="name CA") #if you want choose another atom, just change the name, prefer to this link: https://userguide.mdanalysis.org/stable/selections.html
rmsd.run()

# Use results.rmsd instead of rmsd attribute to avoid deprecation warning
df = pd.DataFrame({
    "Time (ps)": rmsd.results.rmsd[:, 1],
    "RMSD (Å)": rmsd.results.rmsd[:, 2]
})
df.to_csv("rmsd_results.csv", index=False)

# Plot
plt.plot(df["Time (ps)"], df["RMSD (Å)"])
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (Å)")
plt.title("Backbone RMSD")
plt.savefig("rmsd_plot.png")

# Show the plot
plt.show()