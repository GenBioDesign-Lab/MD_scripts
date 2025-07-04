# Step 2: Count water molecules in the solvated system
# This script loads the solvated system and counts water molecules

# Input files from solvation step
set inp_psf "solvated/solvated.psf"
set inp_pdb "solvated/solvated.pdb"

# Load solvated structure
mol new $inp_psf
mol addfile $inp_pdb

# Count water molecules
set water_sel [atomselect top "water"]
set num_water_atoms [$water_sel num]

# Each water molecule has 3 atoms (H2O), so divide by 3
set num_water_molecules [expr $num_water_atoms / 3]

puts "Water statistics:"
puts "- Total water atoms: $num_water_atoms"
puts "- Total water molecules: $num_water_molecules"

# Output the count in a format that can be parsed by the bash script
puts "WATER_COUNT: $num_water_molecules"

# Also get system dimensions for Packmol setup
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

# Export just the solute (non-water) for Packmol
set solute [atomselect top "not water"]
$solute writepdb "solvated/solute_only.pdb"
$solute writepsf "solvated/solute_only.psf"

puts "Solute-only files created:"
puts "- solvated/solute_only.psf"
puts "- solvated/solute_only.pdb"

# Clean up
$water_sel delete
$all delete
$solute delete
mol delete top

exit 