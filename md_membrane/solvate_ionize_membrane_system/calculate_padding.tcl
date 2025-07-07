# calculate_padding.tcl - Calculate solvation padding values
# Usage: source calculate_padding.tcl (requires XSC and atom data)

# Initialize global variables
global final_pad_plus_x final_pad_minus_x final_pad_plus_y final_pad_minus_y final_pad_z

# Fixed Z-direction padding for membrane systems
set fixed_z_padding 15

# Initialize padding variables
set final_pad_plus_x 0
set final_pad_minus_x 0
set final_pad_plus_y 0
set final_pad_minus_y 0
set final_pad_z $fixed_z_padding

# Access global variables from other scripts
global use_xsc_box target_Lx target_Ly target_ox target_oy
global min_x max_x min_y max_y current_span_x current_span_y current_center_x current_center_y

if {$use_xsc_box} {
    # Calculate target box boundaries
    set target_max_x [expr {$target_ox + $target_Lx / 2.0}]
    set target_min_x [expr {$target_ox - $target_Lx / 2.0}]
    set target_max_y [expr {$target_oy + $target_Ly / 2.0}]
    set target_min_y [expr {$target_oy - $target_Ly / 2.0}]
    
    # Calculate asymmetric padding
    set pad_plus_x [expr {$target_max_x - $max_x}]
    set pad_minus_x [expr {$min_x - $target_min_x}]
    set pad_plus_y [expr {$target_max_y - $max_y}]
    set pad_minus_y [expr {$min_y - $target_min_y}]
    
    # Set final values
    set final_pad_plus_x $pad_plus_x
    set final_pad_minus_x $pad_minus_x
    set final_pad_plus_y $pad_plus_y
    set final_pad_minus_y $pad_minus_y
    
} else {
    # Use default symmetric padding
    set final_pad_plus_x -5
    set final_pad_minus_x -5
    set final_pad_plus_y -5
    set final_pad_minus_y -5
} 