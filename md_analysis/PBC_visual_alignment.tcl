# PBC repair
package require pbctools
pbc unwrap -sel all
set cntSelText "segname CNT1 or segname CNT2 or segname CNT3"
pbc wrap -all -compound residue -center com -centersel $cntSelText

#Apply CNT alignment
set original_mol [molinfo top]
set nf [molinfo $original_mol get numframes]

# Set up atom selections for alignment
set cnt_ref [atomselect top $cntSelText frame 0]
set cnt_now [atomselect top $cntSelText]
set all_now [atomselect top "all"]

# Error checking: Make sure CNT selection is not empty
set cnt_count [$cnt_ref num]
if {$cnt_count == 0} {
    puts "ERROR: No atoms found for CNT selection: $cntSelText"
    puts "Check your segment names"
    $cnt_ref delete
    $cnt_now delete
    $all_now delete
    return
}

# Start from frame 1 (frame 0 is already the reference)
for {set f 1} {$f < $nf} {incr f} {
    $cnt_now frame $f
    $cnt_now update
    
    # Calculate transformation to align CNT to reference
    set M [measure fit $cnt_now $cnt_ref]
    
    # Apply transformation to all atoms in this frame
    $all_now frame $f
    $all_now update
    $all_now move $M
}

# Clean up atom selections
$cnt_ref delete
$cnt_now delete
$all_now delete