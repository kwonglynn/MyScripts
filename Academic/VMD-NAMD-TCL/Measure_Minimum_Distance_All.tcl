mol new MIT1_P1_POPA_box.gro
mol addfile PAP1_MD_noWater_noTranslation.xtc waitfor all

set file [open "PAP1_minimum_distances_all.dat" w] 

puts $file "#Time\tMin_Distance"

set N_frames [molinfo top get numframes]

puts $N_frames

for {set f 0} {$f < $N_frames} {incr f} {
        set t [expr {$f/100}]
	puts -nonewline $file "[format %-6d $t]\t"
	set protein [atomselect top "protein" frame $f]	
	set lipid [atomselect top "resname POPA and within 20 of protein" frame $f]
        set dists ""	
	foreach atom_protein [$protein get index] {		
		foreach atom_lipid [$lipid get index] {
			lappend dists [measure bond "$atom_protein $atom_lipid" frame $f]
                }
	}
	set dists [lsort -real $dists]
	set min_dist [lindex $dists 0]
	puts $file "[format %6.2f $min_dist]"
}

flush $file
close $file

exit 0
