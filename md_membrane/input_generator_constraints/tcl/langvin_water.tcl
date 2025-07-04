mol new cnt_mem_ionized.psf
mol addfile cnt_mem_ionized.pdb

set all [atomselect top all]
$all set beta 0.0

set wat [atomselect top "water"]

$wat set beta 1.0

$all writepdb langvin_water.pdb

quit