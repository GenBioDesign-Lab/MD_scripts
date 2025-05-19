#Usage vmd -e z_align.tcl <pdb_file>
#package not installed on vmd 2.0.4 yet, must use on 1.9.3

package require Orient
namespace import Orient::orient

# Try to get input file name from command-line arguments or from the loaded molecule
if {[llength $argv] > 0} {
    set input_pdb [lindex $argv 0]
} else {
    set input_pdb [molinfo top get filename]
}
# Extract base name without extension
set input_base [file rootname [file tail $input_pdb]]
set output_pdb "${input_base}_aligned.pdb"

# Load the input PDB file if not already loaded (optional, see note below)
# mol new $input_pdb

set sel [atomselect top "all"]
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $sel]
$sel writepdb $output_pdb

exit