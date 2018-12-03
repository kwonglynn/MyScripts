import MDAnalysis as mda

u = mda.Universe('complex.prmtop', 'complex.inpcrd')

close_resid = u.select_atoms('name CA and same residue as (around 10 resname LIG)')

residues = list(close_resid.resids)

decom_resid = ','.join(map(str,residues))
decom_resid = '1,' + decom_resid	#Add the resid of the ligand to the list, which is needed by the algorithm

fo = open('06_MMPBSA.in', 'w')
with open('06_MMPBSA_raw.in', 'r') as fi:
    for line in fi:
        if 'XXX' in line:
            line = line.replace('XXX', decom_resid)
        fo.write(line)
