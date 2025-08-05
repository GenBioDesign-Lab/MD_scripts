# Membrane preparation solvation script
# Input file paths
set xsc_file "{xsc_file}"
mol new {system_psf}
mol addfile {system_pdb}

# Create output directories
set tmp_dir "tmp"
set final_dir "final"

if {{![file exists "tmp"]}} {{
    file mkdir "tmp"
}}
if {{![file exists "final"]}} {{
    file mkdir "final"
}}

# Global variables for XSC reading
global target_Lx target_Ly target_Lz target_ox target_oy target_oz use_xsc_box

set target_Lx 0
set target_Ly 0
set target_Lz 0
set target_ox 0
set target_oy 0
set target_oz 0
set use_xsc_box 0

# XSC reading function
proc read_xsc_file {{xsc_file}} {{
    if {{![file exists $xsc_file]}} {{
        return {{}}
    }}
    
    set file_handle [open $xsc_file r]
    set line_count 0
    set box_info {{}}
    
    while {{[gets $file_handle line] >= 0}} {{
        incr line_count
        if {{$line_count == 3}} {{
            set box_info [split $line]
            break
        }}
    }}
    close $file_handle
    
    if {{[llength $box_info] >= 13}} {{
        set target_Lx [lindex $box_info 1]
        set target_Ly [lindex $box_info 5]
        set target_Lz [lindex $box_info 9]
        set target_ox [lindex $box_info 10]
        set target_oy [lindex $box_info 11]
        set target_oz [lindex $box_info 12]
        
        return [list $target_Lx $target_Ly $target_Lz $target_ox $target_oy $target_oz]
    }} else {{
        return {{}}
    }}
}}

# Read XSC if provided
if {{[file exists $xsc_file]}} {{
    set xsc_data [read_xsc_file $xsc_file]
    
    if {{[llength $xsc_data] == 6}} {{
        set target_Lx [lindex $xsc_data 0]
        set target_Ly [lindex $xsc_data 1]
        set target_Lz [lindex $xsc_data 2]
        set target_ox [lindex $xsc_data 3]
        set target_oy [lindex $xsc_data 4]
        set target_oz [lindex $xsc_data 5]
        set use_xsc_box 1
        puts "Using XSC box dimensions: ${{target_Lx}} x ${{target_Ly}} x ${{target_Lz}}"
    }} else {{
        puts "Warning: Could not read XSC file"
        set use_xsc_box 0
    }}
}} else {{
    puts "Warning: XSC file not found, using default padding"
    set use_xsc_box 0
}}

# Load analysis modules
source analyze_atoms.tcl
source calculate_padding.tcl

# Save centered system
set all [atomselect top all]
$all writepsf tmp/system_centered.psf
$all writepdb tmp/system_centered.pdb
$all delete

puts "Padding: X(${{final_pad_minus_x}}, ${{final_pad_plus_x}}) Y(${{final_pad_minus_y}}, ${{final_pad_plus_y}}) Z(±${{final_pad_z}})"

# Solvation and ionization
package require solvate
package require autoionize

solvate tmp/system_centered.psf tmp/system_centered.pdb \\
 -x $final_pad_minus_x -y $final_pad_minus_y -z $final_pad_z \\
 +x $final_pad_plus_x +y $final_pad_plus_y +z $final_pad_z \\
 -o tmp/system_solvated

autoionize -psf tmp/system_solvated.psf -pdb tmp/system_solvated.pdb -sc {salt_concentration} -o final/system_ionized

quit 