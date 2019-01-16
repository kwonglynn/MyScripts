import MDAnalysis as mda
import glob

state_name = glob.glob('state*frame*.pdb')[0]
u = mda.Universe(state_name)

prot_noH = u.select_atoms('protein and not name H*')
prot_noH.write('protein_noH.pdb')

LIG = u.select_atoms('resname LIG')
LIG.write('LIG.pdb')

fo1 = open('complex.pdb', 'w')
fo2 = open('protein-amber.pdb', 'w')
with open('LIG.pdb', 'r') as fi:
    for line in fi:
        if line.startswith('ATOM'):
            fo1.write(line)
    fo1.write("TER\n")

with open('protein_noH.pdb', 'r') as fi:
    for line in fi:
        if line.startswith('ATOM'):
            if 'CD  ILE' in line:
                line = line.replace('CD  ILE', 'CD1 ILE')
            fo1.write(line)
            fo2.write(line)

            if 'CH3 NME' in line:
                fo1.write("TER\n")
                fo2.write("TER\n")
    fo1.write("END\n")
    fo2.write("END\n")

fo1.close()
fo2.close()
