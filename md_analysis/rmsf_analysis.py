import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Get input files either from command line or user input
if len(sys.argv) == 3:
    psf_file = sys.argv[1]
    dcd_file = sys.argv[2]
else:
    print("Usage: python rmsf_analysis.py <path_to_psf_file> <path_to_dcd_file>")
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
ca_atoms = u.select_atoms("name CA")   #if you want choose another atom, just change the name, prefer to this link: https://userguide.mdanalysis.org/stable/selections.html

rmsf = RMSF(ca_atoms)
rmsf.run()

# Use results.rmsf instead of rmsf attribute to avoid deprecation warning
df = pd.DataFrame({
    "Residue": ca_atoms.resids,
    "RMSF (Å)": rmsf.results.rmsf
})
df.to_csv("rmsf_results.csv", index=False)

plt.plot(df["Residue"], df["RMSF (Å)"])
plt.xlabel("Residue")
plt.ylabel("RMSF (Å)")
plt.title("Per-residue RMSF")
plt.savefig("rmsf_plot.png")

# Show the plot
plt.show()