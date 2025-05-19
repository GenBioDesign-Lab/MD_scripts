import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math

# Read the energy.dat file
file_path = r"C:\Users\Blue\Desktop\test\energy.dat"
data = pd.read_csv(file_path, delim_whitespace=True, header=None)

# Print first few rows to see the structure
print("Data preview:")
print(data.head())

# Assuming the energy values are in one of the columns
energy_data = data.iloc[:, 0]

# Create a histogram
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(energy_data, bins=100, range=(0, 10), density=True, alpha=0.7, label='Histogram')

# Calculate bin centers for curve fitting
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width = bins[1] - bins[0]

# Define the Maxwell-Boltzmann-like distribution function
def maxwell_boltzmann(x, a0):
    return (2 / (np.sqrt(np.pi * a0**3))) * np.sqrt(x) * np.exp(-x / a0)

# Curve fitting
# Initial parameter guess
initial_guess = [0.6]  # Initial a0 value (corresponds to ~300K)
try:
    # Only use positive bin values for fitting (avoid zeros that can cause issues)
    mask = (bin_centers > 0) & (counts > 0)
    popt, pcov = curve_fit(maxwell_boltzmann, bin_centers[mask], counts[mask], p0=initial_guess)
    
    # Get the optimized parameter
    a0_fit = popt[0]
    
    # Calculate the correlation coefficient (R^2)
    residuals = counts[mask] - maxwell_boltzmann(bin_centers[mask], a0_fit)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((counts[mask] - np.mean(counts[mask]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Calculate temperature from a0
    kb = 0.0019857  # Boltzmann constant in kcal/(mol*K)
    temperature = a0_fit / kb
    
    # Generate fitted curve for plotting
    x_fit = np.linspace(0.01, 10, 1000)
    y_fit = maxwell_boltzmann(x_fit, a0_fit)
    
    # Plot the fitted curve
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, label=f'Fit: a0 = {a0_fit:.6f}')
        
except Exception as e:
    print(f"Curve fitting error: {e}")

# Add labels and title
plt.xlabel('Energy (kcal/mol)')
plt.ylabel('Normalized Frequency')
plt.title('Maxwell-Boltzmann Energy Distribution')
plt.grid(True, alpha=0.3)
plt.legend()

# Resize to fit
plt.tight_layout()

# Save figure
output_path = r"C:\Users\Blue\Desktop\energy_histogram_fitted.png"
plt.savefig(output_path)
print(f"Figure saved to: {output_path}")

# Show the plot
plt.show()