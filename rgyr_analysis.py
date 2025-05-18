import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Get input files either from command line or user input
if len(sys.argv) == 3:
    psf_file = sys.argv[1]
    dcd_file = sys.argv[2]
else:
    print("Usage: python rgyr_analysis.py <path_to_psf_file> <path_to_dcd_file>")
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
df.to_csv("rgyr_results.csv", index=False)

plt.plot(times, rgyr)
plt.xlabel("Time (ps)")
plt.ylabel("Radius of Gyration (Å)")
plt.title("Protein Compactness")
plt.savefig("rgyr_plot.png")

# Show the plot
plt.show()