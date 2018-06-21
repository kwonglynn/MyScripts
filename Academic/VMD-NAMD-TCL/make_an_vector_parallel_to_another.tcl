# This is the example for pulling DCVJ from alpha-synuclein in the core site.
# Firt put two structures with DCVJ at two different sites in the core site, the ligands
# of which are used to define the vector of pulling.
# References:
# 1. http://www.ks.uiuc.edu/Research/vmd/vmd-1.7/ug/node161.html
# 2. http://www.ks.uiuc.edu/Research/vmd/vmd-1.7/ug/node164.html
# 3. http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/3898.html
# 4. 

mol new state5.gro
mol new state1.gro

set lig1 [atomselect 0 "resname LIG"]
set lig2 [atomselect 1 "resname LIG"]

set COM_lig1 [measure center $lig1]
set COM_lig2 [measure center $lig2]

set vec_lig [vecsub $COM_lig1 $COM_lig2]
set len_vec_lig [veclength $vec_lig]
# The reference axis
set ref_axis {0 0 1}

set prod [vecdot $vec_lig $ref_axis]
# The length for the reference axis vector is unit.
set angle [expr 57.2958 * acos($prod / $len_vec_lig)]

set norm_vec [veccross $vec_lig $ref_axis]

set mat [transabout $norm_vec $angle]

set all [atomselect 0 all]
$all move $mat

##Additional move to make the work easier
#Let the pull route be in the +z direction
set mat2 [trasnaxis y 180]
$all move $mat2

$all writepdb state5_move.pdb
