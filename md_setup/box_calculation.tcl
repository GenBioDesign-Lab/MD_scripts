#Usage: vmd -dispdev text -e box_calculation.tcl <.psf file> <.pdb file>
# Get command-line arguments
set psf_file [lindex $argv 0]
set pdb_file [lindex $argv 1]

# Load system
mol new $psf_file
mol addfile $pdb_file

# Bounding box
set all   [atomselect top all]
set mm    [measure minmax $all]     ;# {{xmin ymin zmin} {xmax ymax zmax}}
set min   [lindex $mm 0]
set max   [lindex $mm 1]

# Box lengths
set Lx [expr {[lindex $max 0] - [lindex $min 0]}]
set Ly [expr {[lindex $max 1] - [lindex $min 1]}]
set Lz [expr {[lindex $max 2] - [lindex $min 2]}]

# Cell origin = midpoint
set origin [vecadd $min [vecscale 0.5 [vecsub $max $min]]]

puts "cell_basis_vector1:  \"$Lx 0.0 0.0\""
puts "cell_basis_vector2:  \"0.0 $Ly 0.0\""
puts "cell_basis_vector3:  \"0.0 0.0 $Lz\""
puts "cell_origin:         \"[format %.3f [lindex $origin 0]] \
                            [format %.3f [lindex $origin 1]] \
                            [format %.3f [lindex $origin 2]]\""
quit
