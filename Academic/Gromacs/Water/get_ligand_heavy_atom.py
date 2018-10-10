# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 18:47:46 2018

@author: Guanglin Kuang
"""

import MDAnalysis as mda

u= mda.Universe('ASEM_move_to_MD.pdb')
noH = u.select_atoms("not name H*")

xyz_all_noH = []
for atom in noH.atoms:
    xyz = atom.position/10 #Change Angstrom to nm
    xyz = list(map(lambda x:format(x, '.2f'), xyz)) #Format the float to have two decimal points
    xyz = ' '.join(xyz)
    xyz_all_noH.append(xyz)

print(','.join(xyz_all_noH))
with open('Ligand_heavy_atoms.xyz', 'w') as fo:
    fo.write(','.join(xyz_all_noH))
    fo.write('\n')
    
##Command line to use Xianqiang's script:
# python job_submission.py -c protein.top -t ../../md_noPBC.xtc -w T3P -gro md.gro --ignore_ww_orientational no --output_pdb water.pdb --volume_pdb volume_water.pdb --output_pickle water.pickle coordinate --center