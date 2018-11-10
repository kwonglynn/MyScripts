import argparse
import MDAnalysis as mda

parser = argparse.ArgumentParser()
parser.add_argument("input", help="The pdb file of the protein.")
parser.add_argument("-o", "--output", default="protein_rename.pdb", help="The name of the output pdb file.")

try:
    args = parser.parse_args()
except:
    print ("\nExample: python path_rename.py protein.pdb -o protein_rename.pdb\n")
    quit()

fi = open(args.input,'r')
fo = open(args.output,'w')

#The atom indices should be the same as those in the gro file for MD simulation.
u = mda.Universe('5KXI-CA-MD.gro')
atom_ids = u.atoms.ids
atom_names = u.atoms.names
res_names = u.atoms.resnames
res_ids = u.atoms.resids
ALIGNs = list(range(185,195)) + list(range(205,211)) 
MSDs = list(range(195,205))

i = 0
for line in fi:
    if line.startswith('HETATM') or line.startswith('ATOM'):
###Common for all PDBs        
        record = line[0:6].strip()
        serial = atom_ids[i]
        name = atom_names[i]
        altLoc = line[16].strip()
        resName = res_names[i]
        chainID = line[21].strip()
        resSeq = res_ids[i]
        iCode = line[26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occu = float(line[54:60].strip())
        temp = float(line[60:66].strip())
        segID = line[72:76].strip()
        element = line[76:78].strip()
        charge = line[78:80].strip()
        i += 1
        ##Standardize atom names        
        if len(name) == 1:
                name = name + ' ' * 2
        elif len(name) == 2:
                name = name + ' ' * 1
        # Generally hydrogens are not used to define path.
        
        #Standardize iCode and resSeq
        if len(iCode.strip()) > 0:
            iCode = ' '
            if len(resSeq.strip()) == 2:
                resSeq = '  ' + str(int(resSeq)+1)
            elif len(resSeq.strip()) == 3:
                resSeq = ' ' + str(int(resSeq)+1)
            elif len(resSeq.strip()) == 4:
                resSeq = str(int(resSeq)+1)

        ## Define the atoms used for alighment and MSD calculation, 
        ## using the occupancy and temperature factor columns, respectively.
        ## For occupancy
        if resSeq in ALIGNs:
            occu = 1.00
        else:
            occu = 0.00
        ## For beta factor
        if resSeq in MSDs:
            temp = 1.00
        else:
            temp = 0.00
        
        fo.write("%-6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"  % \
                 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge)) 
    
fo.write("END\n")
fi.close()
fo.close()
