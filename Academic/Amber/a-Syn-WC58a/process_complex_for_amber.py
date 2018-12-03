import MDAnalysis as mda
import glob

state_name = glob.glob('state*frame*.pdb')[0]
u = mda.Universe(state_name)

prot_noH = u.select_atoms('protein and not name H*')
prot_noH.write('protein_GMX_noH.pdb')

LIG = u.select_atoms('resname LIG')
LIG.write('LIG.pdb')

fo_p = open('protein.pdb', 'w')
fo_c = open('complex.pdb', 'w')

with open('LIG.pdb', 'r') as fi:
    for line in fi:
        if line.startswith('ATOM'):
            fo_c.write(line)
fo_c.write("TER\n")

with open('protein_GMX_noH.pdb', 'r') as fi:
    for line in fi:
        if line.startswith('ATOM'):
            if 'CD  ILE' in line:
                line = line.replace('CD  ILE', 'CD1 ILE')
            fo_p.write(line)
            fo_c.write(line)

            if 'CH3 NME' in line:
                fo_p.write("TER\n")
                fo_c.write("TER\n")

fo_p.write("END\n")
fo_c.write("END\n")

fo_p.close()
fo_c.close()
