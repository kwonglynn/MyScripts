#1, Measure x, y, z dimensions of a box, which will be used to define the x, y and z dimensions of the water box..
set prot [atomselect top protein]
set minmax [measure minmax $prot]
set min [lindex $minmax 0]
set max [lindex $minmax 1]
set dim [vecsub $min $max]

#2， Load a series of coordinate files
for {set x 1} {$x < 71} {incr x} {
    if {[file exists complex_$x.gro]} {
    mol new complex_$x.gro}}
    


