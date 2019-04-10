import pandas as pd
import subprocess

df = pd.read_csv('score_local_dock.sc', skiprows=1, header=0, sep='\s+')

score_min = df['total_score'].min()
id_min = df['total_score'].idxmin()
pdb_min = df.iloc[id_min]['description']

print("Minimum total score: \t{:.2f}".format(score_min))
print("The best PDB structure is: \t{}.pdb".format(pdb_min))

subprocess.call(['cp', pdb_min+'.pdb', '../'+'1000_'+pdb_min+'.pdb'])
