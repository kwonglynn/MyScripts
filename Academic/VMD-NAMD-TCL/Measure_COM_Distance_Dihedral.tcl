#Load and align the strucutre

set file [open "COM_dist_dihed.dat" w]

puts $file "#Frame\tdist_com\tdist\tdihed"

set N_frames [molinfo top get numframes]

for {set f 1} {$f < $N_frames} {incr f} {
        puts -nonewline $file [format "%-6d\t" $f]
        set prot [atomselect top "serial 130 to 146 1590 to 1608" frame $f]
        set protCOM [measure center $prot]
        set pip [atomselect top "serial 1952 to 1978" frame $f]
        set pipCOM [measure center $pip]

        set dist_com [veclength [vecsub $protCOM $pipCOM]]

        set A1 [atomselect top "index 0" frame $f]
        $A1 set x [lindex $pipCOM 0]
        $A1 set y [lindex $pipCOM 1]
        $A1 set z [lindex $pipCOM 2]

        set A2 [atomselect top "index 1" frame $f]
        $A2 set x [lindex $protCOM 0]
        $A2 set y [lindex $protCOM 1]
        $A2 set z [lindex $protCOM 2]

        set dist [measure bond {0 1} frame $f]
        set dihed [measure dihed {1661 1911 0 1} frame $f]
        puts $file [format "%6.2f\t%6.2f\t%6.2f" $dist_com $dist $dihed]
}

flush $file
close $file
