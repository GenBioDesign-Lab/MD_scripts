#!/bin/bash

# Main script to solvate system, count waters, use Packmol for leaflet placement, and clean up
# Usage: bash solvate_packmol_ionize.sh <water_pdb_file>

set -e  # Exit on any error

# Check if water file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <water_pdb_file>"
    echo "Example: $0 /path/to/TIP3.pdb"
    exit 1
fi

WATER_PDB_FILE="$1"

# Function to check if required programs are available
check_dependencies() {
    echo "Checking dependencies..."
    
    if ! command -v vmd &> /dev/null; then
        echo "Error: VMD is not installed or not in PATH"
        exit 1
    fi
    
    if ! command -v packmol &> /dev/null; then
        echo "Error: Packmol is not installed or not in PATH" 
        exit 1
    fi
    
    if ! command -v python3 &> /dev/null; then
        echo "Error: Python3 is not installed or not in PATH"
        exit 1
    fi
    
    echo "All dependencies found."
}

# Configuration
INP_PSF="../cleaned/cnt_mem.psf"
INP_PDB="../cleaned/cnt_mem.pdb"
SOLVATED_DIR="solvated"
PACKMOL_DIR="packmol_output"
IONIZED_DIR="final_ionized"
WATER_THICKNESS=15  # Angstroms for each leaflet

# Check dependencies first
check_dependencies

# Validate input files
if [ ! -f "$INP_PSF" ]; then
    echo "Error: Input PSF file not found: $INP_PSF"
    exit 1
fi

if [ ! -f "$INP_PDB" ]; then
    echo "Error: Input PDB file not found: $INP_PDB"  
    exit 1
fi

if [ ! -f "$WATER_PDB_FILE" ]; then
    echo "Error: Water PDB file not found: $WATER_PDB_FILE"
    exit 1
fi

echo "Input files validated."

# Create directories
mkdir -p $SOLVATED_DIR $PACKMOL_DIR $IONIZED_DIR

echo "=== Step 1: Solvating system with VMD ==="
vmd -dispdev text -e step1_solvate.tcl

echo "=== Step 2: Counting water molecules ==="
WATER_COUNT=$(vmd -dispdev text -e step2_count_water.tcl | grep "WATER_COUNT:" | cut -d: -f2 | tr -d ' ')
echo "Total water molecules found: $WATER_COUNT"

# Calculate waters per leaflet (divide by 2)
WATERS_PER_LEAFLET=$((WATER_COUNT / 2))
echo "Waters per leaflet: $WATERS_PER_LEAFLET"

echo "=== Step 3: Generating Packmol input file ==="
python3 step3_generate_packmol.py $WATERS_PER_LEAFLET $WATER_THICKNESS "$WATER_PDB_FILE"

echo "=== Step 4: Running Packmol ==="
bash step4_run_packmol.sh

echo "=== Cleaning up intermediate files ==="
# Remove solvated directory and files
rm -rf $SOLVATED_DIR
echo "Removed: $SOLVATED_DIR/"

# Remove water template file if it exists in packmol_output
rm -f $PACKMOL_DIR/water.pdb
echo "Removed: $PACKMOL_DIR/water.pdb (if existed)"

# Remove ionized directory (empty anyway since we skipped ionization)
rm -rf $IONIZED_DIR
echo "Removed: $IONIZED_DIR/"

# Keep packmol_leaflets.inp for reproducibility
echo "Kept: packmol_leaflets.inp (Packmol input file)"

echo "=== All steps completed successfully! ==="
echo "Final outputs:"
echo "- System with water leaflets: $PACKMOL_DIR/system_with_leaflets.pdb"
echo "- Packmol input file: packmol_leaflets.inp"
echo "All intermediate files have been cleaned up." 