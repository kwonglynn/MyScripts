#1, Measure x, y, z dimensions of a box, which will be used to define the x, y and z dimensions of the water box..
set prot [atomselect top protein]
set minmax [measure minmax $prot]
set min [lindex $minmax 0]
set max [lindex $minmax 1]
set dim [vecsub $min $max]

