# -*- coding: utf-8 -*-

fi = open('ligand_raw.pdb','r')
fo = open('LIG.pdb','w')

resName =  'LIG' 		# The residue name of the ligand.
rename_atom = 'True'		# Whether to rename the atoms

element2 = ['Cl','CL','Br','BR','PB']

N = 0
for line in fi:
    if 'HETATM' in line or 'ATOM' in line:
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
###ligand specific        
        N += 1
        serial = N

        if rename_atom == 'False':
            if len(name) > 1 and name[:2] in element2:
                if len(name) == 2:
                    name = name + ' ' * 2
                elif len(name) == 3:
                    name = name + ' ' * 1
            else:
                if len(name) == 1:
                    name = name + ' ' * 2
                elif len(name) == 2:
                    name = name + ' ' * 1
        elif rename_atom == 'True':
            if len(name) > 1 and name[:2] in element2:
                name = name[:2] + str(N)
                if len(name) == 3:
                    name = name + ' ' * 1
            else:
                name = name[0] + str(N)
                if len(name) == 2:
                    name = name + ' ' * 1

        resName = 'LIG'
        chainID = 'A'
        resSeq = 1

        fo.write("%-6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"  % \
                 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge))

fo.write("TER\n")

fi.close()
fo.close()
