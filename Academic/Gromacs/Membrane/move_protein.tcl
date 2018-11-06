mol new protein_princ.gro
mol new POPG-noWater-90ns.gro

set prot [atomselect 0 all]
set top_P [atomselect 1 "name P and z > 40"]

set center_prot [measure center $prot]
set center_top_P [measure center $top_P]
set vec_move [vecsub $center_top_P $center_prot]

$prot moveby $vec_move
$prot moveby {0 0 11.5}
$prot writepdb protein_moved.pdb

quit
