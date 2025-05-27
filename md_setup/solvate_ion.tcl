# Usage: vmd -e solvate_ion.tcl -args <input.psf> <input.pdb> <output_dir>
# Parse command-line arguments
set inp_psf [lindex $argv 0]
set inp_pdb [lindex $argv 1]

# Set output directory (default or user-specified)
puts "Number of arguments: $argc"
puts "Arguments: $argv"

if {$argc >= 3} {
    set out_dir [lindex $argv 2]
    puts "Using user-specified output directory: $out_dir"
} else {
    set out_dir "output_ionized"
    puts "Using default output directory: $out_dir"
}

# Extract base name for file prefixes (remove any path components)
set out_prefix [file tail $out_dir]

# Create output directory if it does not exist
file mkdir $out_dir

# Load required plugins
package require solvate
package require autoionize

# 1. Load starting structure
#mol new $inp_psf
#mol addfile $inp_pdb

# 2. SOLVATE
#  -t 10   : add 10 Å padding in every direction
#  -s 1    : use a 1 Å solvent grid spacing (default 1.0 Å; larger = faster)
#  -o ...  : prefix for output files
solvate $inp_psf $inp_pdb -t 10 -s 2 -o $out_dir/${out_prefix}_solv

# The call creates: <out_dir>/<out_prefix>_solv.psf / <out_dir>/<out_prefix>_solv.pdb
# Keep only the solvated molecule in memory
mol delete top
mol load psf $out_dir/${out_prefix}_solv.psf pdb $out_dir/${out_prefix}_solv.pdb

# ----‑‑‑ 3. AUTOIONIZE ------------------------------------------------------
#  -sc       : neutralize and set salt concentration (mol/L)
#  -cation 0.15        : add additional Na⁺/Cl⁻ to reach 0.15 M
#  -o ...            : prefix for ionized output
autoionize -psf $out_dir/${out_prefix}_solv.psf \
           -pdb $out_dir/${out_prefix}_solv.pdb \
           -sc 0.15 \
           -cation SOD \
           -o $out_dir/${out_prefix}_ionized

# Creates: <out_dir>/<out_prefix>_ionized.psf / <out_dir>/<out_prefix>_ionized.pdb

exit