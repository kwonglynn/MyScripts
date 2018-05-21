#!/opt/anaconda2/bin/python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", help="The pdb file of the ligand.")
parser.add_argument("-s", "--start", default=7, type=int, help="The start position of the X coordinate.")
parser.add_argument("-o", "--output", default="LIG.pdb", help="The name of the output pdb file.")

try:
    args = parser.parse_args()
except:
        print '''
Example: clean_LIG.py ligand.pdb -s 6 -o LIG.pdb
'''
        quit()

fi=open(args.input,'r')
fo=open(args.output,'w')
SX=args.start-1
N=0
rename_atom = 'True'
stand_names = ['H','C','O','F','P','S','Cl','CL','Br','BR','I']

for line in fi:
	if 'HETATM' in line or 'ATOM' in line:
		N += 1
		terms=line.strip().split()
		record='ATOM'
		serial=N
		name=terms[2]
		if rename_atom == 'True':
		    if len(name) > 1 and name[:2] in stand_names:
			name = name[:2] + str(N)
		    else:
	                name=name[0]+str(N)
		resName='LIG'
		chain='A'
		resSeq=1
		# Might differ depending on the input file
		x=float(terms[SX])
		y=float(terms[SX+1])
		z=float(terms[SX+2])
#		occu=float(terms[SX+3])
#		temp=float(terms[SX+4])
		occu=float(line[56:60].strip())
		temp=float(line[60:66])
                element=terms[-1]
		fo.write("%-6s%5d %4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"  % \
			 (record, serial, name, resName, chain, resSeq, x, y, z, occu, temp, element))

fo.write("TER\n")

fi.close()
fo.close()
