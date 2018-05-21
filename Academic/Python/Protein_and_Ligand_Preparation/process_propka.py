# -*- coding: utf-8 -*-

fi1=open('a4b2_5KXI_truncated_prep_noH_clean.pqr','r')
fi2=open('a4b2_5KXI_truncated_prep_noH_clean.pdb','r')
fo1=open('Histidines_propka.dat','w')
fo2=open('a4b2_5KXI_truncated_prep_noH_clean_propka.pdb','w')

###Process the pqr file produced by PROPKA.
histidines = []
for line in fi1:
    if 'ATOM' in line:
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
        
        if 'HI' in resName:
            histidine = [resName,chainID,resSeq]
            if histidine not in histidines:
                histidines.append(histidine)

for term in histidines:
    line = term[0]+'\t'+term[1]+'\t'+term[2]+'\n'       
    fo1.write(line)

###Change the histidine lines according to the pqr file.
for line in fi2:
    if 'ATOM' in line:
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
        
        if 'HI' in resName:
            for histidine in histidines:
                if [chainID,resSeq] == histidine[1:]:
                    resName,chainID,resSeq = histidine
                    break
                
        #Process the NME cap residue
        if resName == 'NMA':
            resName = 'NME'
            
            if len(iCode.strip()) > 0:
                iCode = ' '
                resSeq = ' ' + str(int(resSeq)+1)
                
            if name.strip() == 'CA':
                name = ' CH3'
        
        fo2.write("%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"  % \
			 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge))
    elif 'TER' in line or 'END' in line:
        fo2.write(line)    

fi1.close()            
fi2.close()
fo1.close()
fo2.close()
