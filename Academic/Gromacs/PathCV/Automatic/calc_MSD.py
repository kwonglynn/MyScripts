import glob
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from operator import itemgetter

PDBs = glob.glob("*_renum.pdb")

align_region = "name CA and (resid 185:195 or resid 204:210)"
MSD_region = "name CA and (resid 196:203)"

def calc_MSD(ref, mobile, MSD_region):
    xyz_ref = ref.select_atoms(MSD_region).positions
    xyz_mobile = mobile.select_atoms(MSD_region).positions
    MSD = np.square(rms.rmsd(xyz_mobile, xyz_ref))
    return MSD

results = []
ref = mda.Universe('5KXI-A_renum.pdb')
ref_beta_CA = ref.select_atoms(align_region)
N_CA_ref = ref.select_atoms("name CA").n_atoms

for PDB in PDBs:
    code = PDB.split('-')[0]
    mobile = mda.Universe(PDB)

    N_CA_mobile = mobile.select_atoms("name CA").n_atoms
    if N_CA_mobile != N_CA_ref:    continue

    mobile_beta_CA = mobile.select_atoms(align_region)
    RMSD_beta = rms.rmsd(mobile_beta_CA.positions, ref_beta_CA.positions)
    #RMSDs = align.alignto(mobile, ref, select=align_region, weights="mass")
    #RMSD_beta = RMSDs[1]
    MSD_loop = calc_MSD(ref, mobile, MSD_region)
    results.append([code, RMSD_beta, MSD_loop])

results_sorted = sorted(results, key=itemgetter(2))
with open("align_MSD.txt", 'w') as fo:
    fo.write("#Code\tRMSD_beta\tMSD_loopC\n")
    for result in results_sorted:
        fo.write("%s\t%5.2f\t%5.2f\n" % tuple(result))
