import glob
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from operator import itemgetter

#PDBs = glob.glob("*_renum.pdb")
PDBs = {1:'5KXI-A_renum.pdb', 2:'1uv6-5-A_renum.pdb', 3:'5O8T-A_renum.pdb', 4:'2WN9-A_renum.pdb', 5:'2BYN-A_renum.pdb',6:'4BFQ-A_renum.pdb'} 

align_region = "name CA and (resid 185:195 or resid 204:210)"
MSD_region = "name CA and (resid 196:203)"

results = []
for i in range(1,len(PDBs.keys())):
    ref_PDB = PDBs[i]
    mobile_PDB = PDBs[i+1]
    ref = mda.Universe(ref_PDB)
    mobile = mda.Universe(mobile_PDB)

    ref_beta_CA = ref.select_atoms(align_region)
    #mobile_beta_CA = mobile.select_atoms(align_region)
    #RMSD_beta = rms.rmsd(mobile_beta_CA.positions, ref_beta_CA.positions)
    RMSDs = align.alignto(mobile, ref, select=align_region, weights="mass")
    RMSD_beta = RMSDs[1]

    ref_loop_CA = ref.select_atoms(MSD_region)
    mobile_loop_CA = mobile.select_atoms(MSD_region)
    MSD_loop = np.square(rms.rmsd(mobile_loop_CA.positions, ref_loop_CA.positions))
    
    code=str(i) + '-' + str(i+1)
    results.append([code, RMSD_beta, MSD_loop])

#results_sorted = sorted(results, key=itemgetter(2))
with open("align_MSD.txt", 'w') as fo:
    fo.write("#Code\tRMSD_beta\tMSD_loopC\n")
    for result in results:
        fo.write("%s\t%5.2f\t%5.2f\n" % tuple(result))
