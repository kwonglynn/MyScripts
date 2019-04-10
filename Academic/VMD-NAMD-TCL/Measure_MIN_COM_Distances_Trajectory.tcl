#Usage:
#vmd -dispdev text -e Measure_MIN_COM_Distances_Trajectory.tcl

#mol new MIT_POPA_far_box.gro
#mol addfile md1_2.simple.xtc waitfor all

set fo [open "POPA1_MIN_COM_RMSD_DIHED_test.dat" w] 

set reference [atomselect top "backbone" frame 0]

puts $fo "#Time(ns) MIN(nm) COM(nm) RMSD(nm) DIHED1(rad) DIHED2(rad)"

set N_frames [molinfo top get numframes]
for {set f 0} {$f < $N_frames} {set f [expr {$f+20}]} {
## Print the time
    set t [ expr {$f*0.01} ]
    puts -nonewline $fo [format "%-6.2f\t " $t]

## Measure the minimum distance. Fist choose a small section of POPA heavy atoms which are closet to the protein.
    ## Use the heavy atoms!
    set PROTs [atomselect top "noh and protein" frame $f]

    set sep 1
    set flag "true"

    while {$flag=="true"} {
        set POPAs [atomselect top "(noh and resname POPA) and within $sep of (noh and protein)" frame $f ]
        set N_atoms [$POPAs num]
        if {$N_atoms >= 1} {
            set flag "false"
        }  else {
            incr sep
           }
    }

    set dists ""	
    foreach PROT [$PROTs get index] {		
        foreach POPA [$POPAs get index] {
            set dist [measure bond "$PROT $POPA" frame $f] 
            lappend dists $dist
    }
        }

    set dists [lsort -real $dists]
    set dist_MIN [lindex $dists 0]
    set dist_MIN [expr {$dist_MIN*0.1}]
    puts -nonewline $fo [format "%6.2f\t " $dist_MIN]

#Measure COM Distance. Note: the variable names are redefined!
    set PROT [atomselect top "noh and protein" frame $f]
    set COM_PROT [measure center $PROT]
    set POPA [atomselect top "noh and resname POPA" frame $f]
    set COM_POPA [measure center $POPA]
    set dist_COM [veclength [vecsub $COM_PROT $COM_POPA]]
    set dist_COM [expr {$dist_COM*0.1}]
    puts -nonewline $fo [format "%6.2f\t " $dist_COM]

#Measure the RMSD of the protein after alignment of the protein itself
    set compare [atomselect top "backbone" frame $f]
  
    # compute the transformation
    set M [measure fit $compare $reference]
      # do the alignment
    set all [atomselect top "all" frame $f]
    $all move $M

    # compute the RMSD
    set rmsd [measure rmsd $compare $reference]
    set rmsd [expr {$rmsd*0.1}]
    # print the RMSD
    puts -nonewline $fo [format "%6.2f\t  " $rmsd]

#Measure dihedral. Note: the variable names are redefined!
    set L1 [atomselect top "backbone and resid 118 to 121" frame $f]
    set COM_L1 [measure center $L1]

    set L2 [atomselect top "backbone and resid 145 to 148" frame $f]
    set COM_L2 [measure center $L2]

    set PROT [atomselect top "noh and protein" frame $f]
    set COM_PROT [measure center $PROT]

    set POPA [atomselect top "noh and resname POPA" frame $f]
    set COM_POPA [measure center $POPA]

    set v1 [vecsub $COM_L2 $COM_L1]
    set v2 [vecsub $COM_POPA $COM_PROT]
    set dot [vecdot $v1 $v2]
    set lv1 [veclength $v1]
    set lv2 [veclength $v2]
    set dihed1 [expr {acos($dot/$lv1/$lv2)}]
    puts -nonewline $fo [format "%6.2f\t" $dihed1]  

    set A1 [atomselect top "index 1" frame $f]
    $A1 set x [lindex $COM_L1 0]
    $A1 set y [lindex $COM_L1 1]
    $A1 set z [lindex $COM_L1 2]

    set A2 [atomselect top "index 2" frame $f]
    $A2 set x [lindex $COM_L2 0]
    $A2 set y [lindex $COM_L2 1]
    $A2 set z [lindex $COM_L2 2]

    set A3 [atomselect top "index 3" frame $f]
    $A3 set x [lindex $COM_PROT 0]
    $A3 set y [lindex $COM_PROT 1]
    $A3 set z [lindex $COM_PROT 2]

    set A4 [atomselect top "index 4" frame $f]
    $A4 set x [lindex $COM_POPA 0]
    $A4 set y [lindex $COM_POPA 1]
    $A4 set z [lindex $COM_POPA 2]

    set dihed2 [measure dihed {1 2 3 4} frame $f]
    set dihed2 [expr {($dihed2/180)*3.1415926}]
    puts $fo [format "%6.2f" $dihed2]
}

flush $fo
close $fo

#exit
