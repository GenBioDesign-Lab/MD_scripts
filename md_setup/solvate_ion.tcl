# Usage: vmd -dispdev text -e solvate_ion.tcl -args <input.psf> <input.pdb> [output_prefix]
# Parse command-line arguments
set inp_psf [lindex $argv 0]
set inp_pdb [lindex $argv 1]

# Set output prefix (default or user-specified)
puts "Number of arguments: $argc"
puts "Arguments: $argv"

if {$argc >= 3} {
    set out_prefix [lindex $argv 2]
    puts "Using user-specified output prefix: $out_prefix"
} else {
    set out_prefix "protein"
    puts "Using default output prefix: $out_prefix"
}

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
solvate $inp_psf $inp_pdb -t 10 -s 2 -o ${out_prefix}_solv

# The call creates: <out_prefix>_solv.psf / <out_prefix>_solv.pdb
# Keep only the solvated molecule in memory
mol delete top
mol load psf ${out_prefix}_solv.psf pdb ${out_prefix}_solv.pdb

# ----‑‑‑ 3. AUTOIONIZE ------------------------------------------------------
#  -sc       : neutralize and set salt concentration (mol/L)
#  -cation 0.15        : add additional Na⁺/Cl⁻ to reach 0.15 M
#  -o ...            : prefix for ionized output
autoionize -psf ${out_prefix}_solv.psf \
           -pdb ${out_prefix}_solv.pdb \
           -sc 0.15 \
           -cation SOD \
           -o ${out_prefix}_ionized

# Creates: <out_prefix>_ionized.psf / <out_prefix>_ionized.pdb

quit