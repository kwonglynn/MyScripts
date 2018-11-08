mol new protein_princ.gro
mol new POPS-noPBC-100ns-noWater.gro

set prot [atomselect 0 all]
set top_P [atomselect 1 "name P and z > 40"]

set center_prot [measure center $prot]
set center_top_P [measure center $top_P]
set vec_move [vecsub $center_top_P $center_prot]

$prot moveby $vec_move
$prot moveby {0 0 12}
$prot writepdb protein_P1.pdb

set center0 [measure center $prot]
$prot move [transaxis x 90]
set center1 [measure center $prot]
set vec [vecsub $center0 $center1]
$prot moveby $vec
$prot moveby {0 0 2}
$prot writepdb protein_P2.pdb

set center0 [measure center $prot]
$prot move [transaxis x 90]
set center1 [measure center $prot]
set vec [vecsub $center0 $center1]
$prot moveby $vec
$prot writepdb protein_P3.pdb

set center0 [measure center $prot]
$prot move [transaxis x 90]
set center1 [measure center $prot]
set vec [vecsub $center0 $center1]
$prot moveby $vec
$prot moveby {0 0 -2}
$prot writepdb protein_P4.pdb

quit
