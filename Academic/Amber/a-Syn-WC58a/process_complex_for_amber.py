import MDAnalysis as mda
import glob

state_name = glob.glob('state*frame*.pdb')[0]
u = mda.Universe(state_name)

prot_noH = u.select_atoms('protein and not name H*')
prot_noH.write('protein_noH.pdb')

LIG = u.select_atoms('resname LIG')
LIG.write('LIG.pdb')

fo = open('complex.pdb', 'w')
with open('LIG.pdb', 'r') as fi:
    for line in fi:
        if line.startswith('ATOM'):
            fo.write(line)
    fo.write("TER\n")

with open('protein_noH.pdb', 'r') as fi:
    for line in fi:
        if line.startswith('ATOM'):
            if 'CD  ILE' in line:
                line = line.replace('CD  ILE', 'CD1 ILE')
            fo.write(line)

            if 'CH3 NME' in line:
                fo.write("TER\n")
    fo.write("END\n")

fo.close()
