import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="The pdb file of the protein.")
parser.add_argument("-o", "--output", default="protein_prep_NME.pdb", help="The name of the output pdb file.")

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
        record = line[0:6].strip()
        serial = line[6:11].strip()
        name = line[12:16].strip()
        altLoc = line[16].strip()
        resName = line[17:20].strip()
        chainID = line[21].strip()
        resSeq = line[22:26].strip()
        iCode = line[26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occu = float(line[54:60].strip())
        temp = float(line[60:66].strip())
        segID = line[72:76].strip()
        element = line[76:78].strip()
        charge = line[78:80].strip()
        
#Protein specific
        #Process the NME cap residue
        if resName == 'NMA':
            resName = 'NME'
            
            if len(iCode.strip()) > 0:
                iCode = ' '
                resSeq = ' ' + str(int(resSeq)+1)
                
            if name.strip() == 'CA':
                name = ' CH3'
        
        fo.write("%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"  % \
			 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge))
    elif 'TER' in line or 'END' in line:
        fo.write(line)