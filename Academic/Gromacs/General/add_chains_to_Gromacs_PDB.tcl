puts "Usage:"
puts "vmd -dispdev text -e add_chains_to_Gromacs_PDB.tcl -args input output" 

set input [lindex $argv 0]
puts "Input is $input"

set output [lindex $argv 1]
puts "Output is $output"

mol new $input

set all [atomselect top "all"]
set TOTAL 75
set chains {A B C D E}

set N_residues [llength [lsort -unique [$all get residue]]]

for {set i 1} {$i < $N_residues} {incr i} {
    set c [expr {($i-1) / $TOTAL}]
    set sel [atomselect top [format "residue %d" $i]]
    set res [$sel get residue]
    set chain [lindex $chains $c]
    $sel set chain $chain
}

set LIG [atomselect top "resname LIG"]
$LIG writepdb LIG.pdb

set prot_no_CAP [atomselect top "protein and not resname ACE NME"]
$prot_no_CAP writepdb $output

quit
