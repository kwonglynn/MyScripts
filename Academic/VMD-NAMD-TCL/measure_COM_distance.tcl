#Load and align the strucutre

set file [open "CV1_CV2.dat" w]

puts $file "#Frame\tdist_com\tdist\t%dihed"

set N_frames [molinfo top get numframes]

for {set f 1} {$f < $N_frames} {incr f} {
        puts -nonewline $file [format "%-6d\t" $f]
        set prot [atomselect top "protein and resid 83 168" frame $f]
        set protCOM [measure center $prot]
        set pip [atomselect top "resname PIP" frame $f]
        set pipCOM [measure center $pip]

        set dist_com [veclength [vecsub $protCOM $pipCOM]]

        set A1 [atomselect top "index 1" frame $f]
        $A1 set x [lindex $protCOM 0]
        $A1 set y [lindex $protCOM 1]
        $A1 set z [lindex $protCOM 2]

        set A2 [atomselect top "index 2" frame $f]
        $A2 set x [lindex $pipCOM 0]
        $A2 set y [lindex $pipCOM 1]
        $A2 set z [lindex $pipCOM 2]

        set dist [measure bond {1 2} frame $f]

	set r173 [atomselect top "name CA and resid 173" frame $f]
        set r188 [atomselect top "name CA and resid 188" frame $f]

        set A3 [atomselect top "index 3" frame $f]
        $A3 set x [$r173 get x]
        $A3 set y [$r173 get y]
        $A3 set z [$r173 get z]

        set A4 [atomselect top "index 4" frame $f]
        $A4 set x [$r188 get x]
        $A4 set y [$r188 get y]
        $A4 set z [$r188 get z]

	set dihed [measure dihed {3 4 2 1} frame $f]
	set dihed [ expr { $dihed/180*3.14 } ]

        puts $file [format "%6.2f\t%6.2f\t%6.2f" $dist_com $dist $dihed]
}

flush $file
close $file
