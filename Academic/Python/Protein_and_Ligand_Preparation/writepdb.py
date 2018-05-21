# -*- coding: utf-8 -*-

fi1=open('a4b2_5KXI_truncated_prep_noH_clean.pdb','r')
#fi2='a4b2_5KXI_truncated_prep_noH_clean.pqr'
fo=open('a4b2_5KXI_truncated_prep_noH_clean_propka.pdb','w')

for line in fi1:
    if 'HETATM' in line or 'ATOM' in line:
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
        fo.write("%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"  % \
			 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge))
    elif 'TER' in line or 'END' in line:
        fo.write(line)    

fi1.close()            
fo.close()