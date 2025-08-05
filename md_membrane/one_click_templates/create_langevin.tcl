mol new {system_psf}
mol addfile {system_pdb}

set all [atomselect top all]
$all set beta 0.0

set wat [atomselect top "water"]
$wat set beta 1.0

$all writepdb {output_pdb}

quit 