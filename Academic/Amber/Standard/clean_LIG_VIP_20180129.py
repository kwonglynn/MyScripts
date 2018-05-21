#!/opt/anaconda2/bin/python

# User Input Here:
resName =  'LIG' 		# The residue name of the ligand.
rename_atom = 'True'		# Whether to rename the atoms
########################################################################
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", help="The pdb file of the ligand.")
parser.add_argument("-o", "--output", default="LIG.pdb", help="The name of the output pdb file.")

try:
    args = parser.parse_args()
except:
        print '''
Example: clean_LIG.py ligand.pdb -o LIG.pdb
'''
        quit()

fi = open(args.input,'r')
fo = open(args.output,'w')

element2 = ['Cl','CL','Br','BR','PB']

if rename_atom: N = 0
for line in fi:
	if 'HETATM' in line or 'ATOM' in line:
		N += 1
		record = line[0:7].strip()
		serial = N
		atName = line[12:16].strip()

                if rename_atom == 'False':
                    if len(atName) > 1 and atName[:2] in element2:
			if len(atName) == 2:
			    atName = atName + ' ' * 2
			elif len(atName) == 3:
			    atName = atName + ' ' * 1
                    else:
			if len(atName) == 1:
			    atName = atName + ' ' * 2
			elif len(atName) == 2:
			    atName = atName + ' ' * 1
		elif rename_atom == 'True':
		    if len(atName) > 1 and atName[:2] in element2:
			atName = atName[:2] + str(N)
			if len(atName) == 3:
			   atName = atName + ' ' * 1
		    else:
	                atName = atName[0] + str(N)
			if len(atName) == 2:
			    atName = atName + ' ' * 1

		chain='A'
		resSeq=1
		# Might differ depending on the input file
		x = float(line[30:38].strip())
		y = float(line[38:46].strip())
		z = float(line[46:54].strip())
		occu = float(line[54:60].strip())
		temp = float(line[60:66].strip())
                element = line[76:78]
		charge = line[78:80]
		altLoc = ' ' * 1
		insertCode = ' ' * 1
		segID = ' ' * 4
		fo.write("%-6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"  % \
			 (record,serial,atName,altLoc,resName,chain,resSeq,insertCode,x,y,z,occu,temp,segID,element,charge))

fo.write("TER\n")

fi.close()
fo.close()
