# -*- coding: utf-8 -*-

"""
Author:
    Guanglin Kuang <guanglin@kth.se>

Usage:
    protein_propka.py [options]

Options:
    --pdb <pdb>                 The PDB file of the protein [default: protein_prep.pdb]
    --his <his>                 The PQR file prepared by ProPKA [default: histidines.dat]
    -o, --output <file>         Save the plot to a file [default: protein_prep_HIS_NME.pdb].

"""
from docopt import docopt
opts = docopt(__doc__)

f_pdb = open(opts["--pdb"], 'r')
f_his = open(opts["--his"], 'r')
fo = open(opts["--output"], 'w')

###Process the pqr file produced by PROPKA.
histidines = []
for line in f_his:
    if (not len(line.strip()) == 0) and (not line.startswith('#')):
       histidines.append(line.strip().split())

###Change the histidine lines according to the pqr file.
for line in f_pdb:
    if line.startswith('ATOM'):
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
                if [resName.strip(), chainID.strip(), resSeq.strip()] == histidine:
                    resName = histidine[0]
                    break
                
        #Process the NME cap residue
        if resName == 'NMA':
            resName = 'NME'
            
            if len(iCode.strip()) > 0:
                iCode = ' '
                if len(resSeq.strip()) == 2:
                    resSeq = '  ' + str(int(resSeq)+1)
                elif len(resSeq.strip()) == 3:
                    resSeq = ' ' + str(int(resSeq)+1)
                elif len(resSeq.strip()) == 4:
                    resSeq = str(int(resSeq)+1)
                
            if name.strip() == 'CA':
                name = ' CH3'
        
        fo.write("%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"  % \
			 (record,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occu,temp,segID,element,charge))

    elif line.startswith('TER') or line.startswith('END'):
        if line.startswith('ENDMDL'):
            continue
        elif not 'NMA' in line:
            fo.write(line)
        elif 'NMA' in line:
            #Common part as before.
            record = line[0:6]
            serial = line[6:11]
            name = line[12:16]
            altLoc = line[16]
            resName = line[17:20]
            chainID = line[21]
            resSeq = line[22:26]
            iCode = line[26]
            
            #Specific for NME
            resName = 'NME'            
            if len(iCode.strip()) > 0:
                iCode = ' '
                if len(resSeq.strip()) == 2:
                    resSeq = '  ' + str(int(resSeq)+1)
                elif len(resSeq.strip()) == 3:
                    resSeq = ' ' + str(int(resSeq)+1)
                elif len(resSeq.strip()) == 4:
                    resSeq = str(int(resSeq)+1)

            fo.write("%s%s %s%s%s %s%s%s\n"  % \
    			 (record,serial,name,altLoc,resName,chainID,resSeq,iCode))
                
f_pdb.close()            
f_his.close()
fo.close()
