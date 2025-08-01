#load file
mol new cnt_mem_ionized.psf
mol addfile cnt_mem_ionized.pdb

#set all atom beta =0

set all [atomselect top all]
$all set beta 0.0

#set k value

set k 2.5 ;# Kcal/mol/A^2

#select cnt

set cnt [atomselect top "segname CNT"]

#store k value in B column
$cnt set beta $k

$all writepdb cnt_constraints.pdb

quit