# analyze_atoms.tcl - Analyze atom coordinates and set global variables
# Usage: source analyze_atoms.tcl (requires loaded molecule)

# Initialize global variables
global min_x max_x min_y max_y min_z max_z current_span_x current_span_y current_span_z current_center_x current_center_y current_center_z

# Get all atoms
set all [atomselect top all]

# Calculate current center of mass
set com [measure center $all weight mass]
set com_x [lindex $com 0]
set com_y [lindex $com 1] 
set com_z [lindex $com 2]

# Move system to origin (center at 0,0,0)
set move_vector [list [expr {-$com_x}] [expr {-$com_y}] [expr {-$com_z}]]
$all moveby $move_vector

# Get updated coordinates after centering
set coords [$all get {x y z}]

# Initialize min/max values
set first_coord [lindex $coords 0]
set min_x [lindex $first_coord 0]
set max_x [lindex $first_coord 0]
set min_y [lindex $first_coord 1]
set max_y [lindex $first_coord 1]
set min_z [lindex $first_coord 2]
set max_z [lindex $first_coord 2]

# Find actual min/max coordinates after centering
foreach coord $coords {
    set x [lindex $coord 0]
    set y [lindex $coord 1]
    set z [lindex $coord 2]
    
    if {$x < $min_x} {set min_x $x}
    if {$x > $max_x} {set max_x $x}
    if {$y < $min_y} {set min_y $y}
    if {$y > $max_y} {set max_y $y}
    if {$z < $min_z} {set min_z $z}
    if {$z > $max_z} {set max_z $z}
}

# Calculate system span and center (should be ~0,0,0 now)
set current_span_x [expr {$max_x - $min_x}]
set current_span_y [expr {$max_y - $min_y}]
set current_span_z [expr {$max_z - $min_z}]

set current_center_x [expr {($min_x + $max_x) / 2.0}]
set current_center_y [expr {($min_y + $max_y) / 2.0}]
set current_center_z [expr {($min_z + $max_z) / 2.0}]

# Clean up
$all delete 