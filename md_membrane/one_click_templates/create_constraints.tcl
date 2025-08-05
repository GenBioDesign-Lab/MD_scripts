mol new {system_psf}
mol addfile {system_pdb}

set all [atomselect top all]
$all set beta 0.0

set k {constraint_force}
set cnt [atomselect top "segname CNT"]
$cnt set beta $k

$all writepdb {output_pdb}

quit 