# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 19:09:30 2018

@author: Guanglin Kuang
"""

mol new md.gro
mol new protein.gro
mol new ASEM_State5.pdb
set sel0 [atomselect 0 "name CA and serial 1 to 3298"]
set sel1 [atomselect 1 "name CA and serial 1 to 3298"]
set sel2 [atomselect 2 all]
set M [measure fit $sel1 $sel0]
$sel2 move $M
$sel2 writepdb Ligand_align_move_to_apo.pdb
quit
