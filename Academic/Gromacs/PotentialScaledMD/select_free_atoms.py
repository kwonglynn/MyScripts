# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np

fo = open('Abeta-Distances_to_COM.dat','w')

u = mda.Universe('npt.gro')

COM = u.select_atoms("bynum 536")

protein = u.select_atoms("protein")

for atom in protein:
   VEC = atom.position - COM.centroid()
   DIST = np.linalg.norm(VEC)
   fo.write("%d\t%4.1f\n" % (atom.index, DIST))

fo.close()