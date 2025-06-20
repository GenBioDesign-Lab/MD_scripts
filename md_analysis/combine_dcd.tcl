#Put this script into the namd folder

for {set i 1} {$i <= 6} {incr i} {
    catdcd -o eq_all.dcd eq$i/final_eq$i.dcd
}

exit