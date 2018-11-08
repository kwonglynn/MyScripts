import numpy as np
import mdtraj as md
from itertools import combinations


NATIVE_CUTOFF = 0.45  # nanometers

traj = md.load_pdb('protein.pdb')

# get the indices of all of the heavy atoms 
heavy = traj.topology.select_atom_indices('heavy')
# get the pairs of heavy atoms which are farther than 3 # residues apart 
heavy_pairs = np.array(
    [(i,j) for (i,j) in combinations(heavy, 2)
        if abs(traj.topology.atom(i).residue.index - \
               traj.topology.atom(j).residue.index) > 3])

# compute the distances between these pairs in the native state 
heavy_pairs_distances = md.compute_distances(traj[0], heavy_pairs)[0] # and get the pairs s.t. the distance is less than NATIVE_CUTOFF 
native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF] # choose th atom number according to chian 
chain_contact = np.array([x for x in native_contacts if x[0] <= 319]) 
chain_contact_distance = md.compute_distances(traj[0], chain_contact)[0]

number_contact = chain_contact.shape[0]
weight = float(1.0) / float(number_contact)

print "CONTACTMAP ..."
for y in range(0,number_contact):
	print("ATOMS%s=%s,%s SWITCH%s={Q R_0=0.01 BETA=50.0 LAMBDA=1.8 REF=%.6f} WEIGHT%s=%.6f" % (y+1, chain_contact[y][0]+1,chain_contact[y][1]+1,y+1,chain_contact_distance[y],y+1,weight))
print "LABEL=cmap"
print "SUM"
print "... CONTACTMAP"

