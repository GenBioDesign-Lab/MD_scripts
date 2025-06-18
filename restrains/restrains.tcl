########## adjustable constants ##########
set outfile       "restraints_POPC.dihe"
set phi0_dihe     0.0
set phi0_glycerol -120.0
set phi0_phosphate 0.0
set tail1         {C21 C22 C23 C24}
set tail2         {C31 C32 C33 C34}
set glycerol_improper {C1 C2 C3 O21}
set phosphate_improper {O11 P O12 O13}
###############################################

set dihedralList [list $tail1 $tail2]
set improperList [list $glycerol_improper $phosphate_improper]

set selAllPOPC  [atomselect top "resname POPC"]
set residList   [lsort -integer -unique [$selAllPOPC get resid]]
$selAllPOPC delete

set fh [open $outfile w]

foreach resid $residList {
    foreach fourAtoms $dihedralList {
        array unset idx
        foreach atomName $fourAtoms {
            set sel [atomselect top "resid $resid and resname POPC and name $atomName"]
            set idx($atomName) [lindex [$sel get index] 0]
            $sel delete
        }
        puts $fh "DIHEDRAL $idx([lindex $fourAtoms 0]) $idx([lindex $fourAtoms 1]) $idx([lindex $fourAtoms 2]) $idx([lindex $fourAtoms 3]) \$FC $phi0_dihe"
    }

    foreach fourAtoms $improperList {
        array unset idx
        foreach atomName $fourAtoms {
            set sel [atomselect top "resid $resid and resname POPC and name $atomName"]
            set idx($atomName) [lindex [$sel get index] 0]
            $sel delete
        }
        if {$fourAtoms eq $glycerol_improper} {
            set phi0 $phi0_glycerol
        } else {
            set phi0 $phi0_phosphate
        }
        puts $fh "IMPROPER $idx([lindex $fourAtoms 0]) $idx([lindex $fourAtoms 1]) $idx([lindex $fourAtoms 2]) $idx([lindex $fourAtoms 3]) \$FC $phi0"
    }
}

close $fh
puts ">> wrote restraints to $outfile" 