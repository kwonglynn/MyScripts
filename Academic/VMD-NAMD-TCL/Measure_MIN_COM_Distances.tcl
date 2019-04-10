set file [open "MIN_COM_distances.dat" w] 

## Use the heavy atoms!
set PROTs [atomselect top "noh and protein"]

puts $file "Min_Distance\tCOM_Distance"

## Measure the minimum distance. Fist choose a small section of POPC heavy atoms which are closet to the protein.
set sep 1
set flag "true"
while {$flag=="true"} {
    set POPCs [atomselect top "(noh and resname POPC) and within $sep of (noh and protein)"]
    set N_atoms [$POPCs num]
    if {$N_atoms >= 1} {
        set flag "false"
        puts "The distance is $sep"
        puts "The number of close atoms is $N_atoms"
    }  else {
          incr sep
    }
}

set dists ""	
foreach PROT [$PROTs get index] {		
    foreach POPC [$POPCs get index] {
        set dist [measure bond "$PROT $POPC"] 
        lappend dists $dist
    }
}

set dists [lsort -real $dists]
set dist_MIN [lindex $dists 0]
puts -nonewline $file [format "%6.2f\t" $dist_MIN]

#Measure COM Distance. Note: the variable names are redefined!
set PROT [atomselect top "protein"]
set COM_PROT [measure center $PROT]
set POPC [atomselect top "resname POPC"]
set COM_POPC [measure center $POPC]
set dist_COM [veclength [vecsub $COM_PROT $COM_POPC]]
puts $file [format %6.2f $dist_COM]

flush $file
close $file
