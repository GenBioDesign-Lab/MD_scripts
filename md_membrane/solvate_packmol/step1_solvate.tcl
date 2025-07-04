# Step 1: Solvate the system using VMD
# This script loads the structure, solvates it, and saves the solvated system

# Input files
set inp_psf "../cleaned/cnt_mem.psf"
set inp_pdb "../cleaned/cnt_mem.pdb"
set out_dir "solvated"

# Load required packages
package require solvate

# Create output directory
file mkdir $out_dir

# Load starting structure
mol new $inp_psf
mol addfile $inp_pdb

# Get system dimensions to determine solvation box
set all [atomselect top "all"]
set minmax [measure minmax $all]
set min_coords [lindex $minmax 0]
set max_coords [lindex $minmax 1]

set min_x [lindex $min_coords 0]
set min_y [lindex $min_coords 1] 
set min_z [lindex $min_coords 2]
set max_x [lindex $max_coords 0]
set max_y [lindex $max_coords 1]
set max_z [lindex $max_coords 2]

puts "System dimensions:"
puts "X: $min_x to $max_x"
puts "Y: $min_y to $max_y" 
puts "Z: $min_z to $max_z"

# Solvate with appropriate padding
# Add 20 Å padding in all directions to ensure good solvation
solvate $inp_psf $inp_pdb \
    -x 0 -y 0  -z 15 \
    +x 0 +y 0  +z 15 \
    -o $out_dir/solvated

puts "Solvation completed. Files saved:"
puts "- $out_dir/solvated.psf"
puts "- $out_dir/solvated.pdb"

# Clean up
$all delete
mol delete top

exit 