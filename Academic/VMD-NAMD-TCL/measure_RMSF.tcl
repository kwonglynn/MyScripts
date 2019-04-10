set f [open "MIT_RMSF_CA_MD1.dat" w]
set atomsCA [atomselect top "name CA"]
set rmsf [measure rmsf $atomsCA]
set IDs [$atomsCA get resid]
set n [llength $rmsf]

for {set i 0} {$i<$n} {incr i} {
    set ID [lindex $IDs $i]
    puts -nonewline $f "$ID\t"
    puts $f [lindex $rmsf $i]
}

close $f
