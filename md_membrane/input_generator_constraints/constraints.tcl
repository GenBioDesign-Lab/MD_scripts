#load file
mol new complex.psf
mol addfile complex.pdb

#set k value

set k 2.5 ;# Kcal/mol/A^2

#select cnt

set cnt [atomselect top "segname CNT"]

#store k value in B column
$cnt set beta $k

$cnt writepdb cnt_constraints.pdb

quit