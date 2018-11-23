import MDAnalysis as mda

u = mda.Universe('complex.prmtop', 'complex.inpcrd')

close_resid = u.select_atoms('name CA and same residue as (around 10 resname LIG)')

residues = list(close_resid.resids)

decom_resid = ','.join(map(str,residues))

fo = open('06_MMPBSA.in', 'w')
with open('06_MMPBSA_raw.in', 'r') as fi:
    for line in fi:
        if 'XXX' in line:
            line = line.replace('XXX', decom_resid)
        fo.write(line)
