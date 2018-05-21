# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import Bio.PDB
import numpy
import pylab

pdblist = Bio.PDB.PDBList()
pdblist.retrieve_pdb_file('1XI4', pdir='.', file_format='pdb')

pdb_code = "1XI4"
pdb_filename = "1XI4.pdb"

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]

dist_matrix = calc_dist_matrix(model["D"], model["M"])
contact_map = dist_matrix < 12.0

pylab.matshow(numpy.transpose(dist_matrix))
pylab.colorbar()
pylab.show()

#pylab.autumn()
#pylab.imshow(numpy.transpose(contact_map))
#pylab.show()