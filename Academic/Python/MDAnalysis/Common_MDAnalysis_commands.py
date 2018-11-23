# -*- coding: utf-8 -*-
from __future__ import print_function
import MDAnalysis as mda

#1. Print the CA atoms around the ligand
u = mda.Universe('complex.gro')
ref = u.select_atoms('(name CA) and (around 8 resname LIG) and (bynum 4499:5392)')
print (*ref.ids)

#2. 
