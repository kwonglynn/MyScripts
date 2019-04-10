###Preprocess of the trajectory with Gromacs
# make_ndx_mpi -f md1.gro -o index.ndx (Combine protein and lipids)
# trjconv_mpi -s md1.tpr -f md1.gro -pbc cluster -ur compact -n index.ndx -o PIP3_Random_MD_A.noPBC.noWater.gro
# trjconv_mpi -s md1.tpr -f md1.xtc -pbc cluster -ur compact -n index.ndx -o PIP3_Random_MD_A.noPBC.noWater.xtc 

set file [open "PIP3_PH_MIN_dist.dat" w]

puts $file "# Time\tMin_Distance"

set N_frames [molinfo top get numframes]

for {set f 1} {$f < $N_frames} {incr f} {
	set time [expr {$f * 0.01}]
	puts -nonewline $file "[format %-7.2f $time]\t"
	set PIP3s [atomselect top "resname PIP3" frame $f]
        
	## Choose a small section of atoms which are closet to the ligand.
	set sep 1
	set flag "true"
        while {$flag=="true"} {
		set PROTs [atomselect top "protein and within $sep of resname PIP3" frame $f]
		set N_atoms [$PROTs num]
 		if {$N_atoms > 1} {
			set flag "false"
			puts "At $time ns"
			puts "The distance is $sep"
			puts "The number of close atoms is $N_atoms"
		} else {
			incr sep
		}
	}
	
        set dists ""	
	foreach PIP3 [$PIP3s get index] {		
		foreach PROT [$PROTs get index] {
			lappend dists [measure bond "$PIP3 $PROT" frame $f]
                }
	}
	set dists [lsort -real $dists]
	puts $dists
	set min_dist [lindex $dists 0]
	puts $min_dist
	puts $file "[format %6.2f $min_dist]"
}

flush $file
close $file
