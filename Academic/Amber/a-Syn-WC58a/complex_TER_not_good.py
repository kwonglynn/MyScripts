import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="The pdb file of the complex.")
parser.add_argument("-o", "--output", default="complex.pdb", help="The name of the output pdb file.")

try:
    args = parser.parse_args()
except:
    print ("\nExample: python protein_NME.py protein_prep.pdb -o protein_prep_NME.pdb\n")
    quit()

fi = open(args.input,'r')
fo = open(args.output,'w')

for line in fi:
    if line.startswith('HETATM') or line.startswith('ATOM'):
###Common for all PDBs        
        record = line[0:6]
        serial = line[6:11]
        name = line[12:16]
        altLoc = line[16]
        resName = line[17:20]
        chainID = line[21]
        resSeq = line[22:26]
        iCode = line[26]
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        occu = line[54:60]
        temp = line[60:66]
        segID = line[72:76]
        element = line[76:78]
        charge = line[78:80]
        
        fo.write("%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s"  % \
			 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge))

        #Add an 'TER' record after the NME residue
        if 'HH33 NME' in line or 'H37 LIG' in line:
            fo.write("TER\n")

    elif 'TER' in line or 'END' in line:
        fo.write(line)

fi.close()
fo.close()
