#!/bin/bash

# Step 4: Run Packmol with error checking
# This script runs Packmol and checks for successful completion

set -e

echo "=== Running Packmol ==="

# Check if Packmol input file exists
if [ ! -f "packmol_leaflets.inp" ]; then
    echo "Error: Packmol input file 'packmol_leaflets.inp' not found!"
    exit 1
fi

# Check if Packmol is available
if ! command -v packmol &> /dev/null; then
    echo "Error: Packmol is not installed or not in PATH"
    echo "Please install Packmol first"
    exit 1
fi

# Create output directory
mkdir -p packmol_output

# Run Packmol
echo "Running Packmol..."
packmol < packmol_leaflets.inp

# Check if output was created
if [ ! -f "packmol_output/system_with_leaflets.pdb" ]; then
    echo "Error: Packmol did not produce expected output file"
    echo "Check packmol_leaflets.inp for errors"
    exit 1
fi

echo "Packmol completed successfully!"
echo "Output file: packmol_output/system_with_leaflets.pdb"

# Get some basic statistics
num_lines=$(wc -l < packmol_output/system_with_leaflets.pdb)
num_atoms=$(grep -c "^ATOM\|^HETATM" packmol_output/system_with_leaflets.pdb || echo 0)

echo "Statistics:"
echo "- Total lines in PDB: $num_lines"
echo "- Total atoms: $num_atoms" 