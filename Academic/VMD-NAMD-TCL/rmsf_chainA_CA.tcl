set f [open "2BEG_RMSF_CA_A.dat" w]
set atomsCA [atomselect top "name CA and chain A"]
set rmsf [measure rmsf $atomsCA]
set IDs [$atomsCA get resid]
set n [llength $rmsf]

for {set i 0} {$i<$n} {incr i} {
    set ID [lindex $IDs $i]
    puts -nonewline $f "$ID\t"
    puts $f [lindex $rmsf $i]
}

close $f
