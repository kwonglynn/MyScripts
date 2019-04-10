set f [open "2BEG_RMSF.dat" w]
set all [atomselect top all]
set rmsf [measure rmsf $all]
set n [llength $rmsf]
for {set i 0} {$i<$n} {incr i} {
puts -nonewline $f "$i\t"
puts $f [lindex $rmsf $i]
}

