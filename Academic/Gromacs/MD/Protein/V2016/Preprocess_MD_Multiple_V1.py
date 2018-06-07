import os
#############################################################################################
#User input
#Protein coordinate file (In)
#Note: Normally only change here!!!
coorP = "2N0A_processed.gro"
#Ligand coordinate file (Check the naming with the topology file, must be consistent!!!) (In)
coorL = "LIG_MD.gro"
#Complex coordinate (Out)
coorC = "complex.gro"

#Ligand name (Molecule type)
nameL = 'LIG'
#Ligand topology file generated by antechamber and acepype (In)
topL0 = "LIG_GMX.top"
#Ligand output topology file for classical MD (Out)
topLMD = "LIG_MD.top"

#The whole topology of protein with multiple chains (In)
top = "topol.top"
#Protein output topology file for MD, the whole one (Out)
topMD = "complex.top"

##########################################################################################
#Combine the coordinate files of protein and ligand.
fcoorP = open(coorP,'r')
fcoorL = open(coorL,'r')
fcoorC = open(coorC,'w')

coorPs = fcoorP.readlines()
coorLs = fcoorL.readlines()

NP = int(coorPs[1].strip())
NL = int(coorLs[1].strip())
NC = NP + NL

fcoorC.write(coorPs[0])
fcoorC.write("%d\n" % NC)

for i in range(2,len(coorLs)-1):
    fcoorC.write(coorLs[i])

for j in range(2,len(coorPs)):
        fcoorC.write(coorPs[j])

fcoorP.close()
fcoorL.close()
fcoorC.close()

#############################################################################################
#Process the ligand topology file.
#The blank lines of the topology files are very important! Don't change the default settings!
#Topology file is generated by antechamber and amb2gmx.pl
ftopL0 = open(topL0,'r')
ftopLMD = open(topLMD,'w')
fatomtypesMD = open('LIG_atomtypes_MD.top','w')

flag = 'NULL'
for line in ftopL0:
    #Process the [ atomtypes ] section.
    if line.startswith('[ atomtypes ]'):
        flag = 'atomtypes'
	fatomtypesMD.write(line)
    elif flag == 'atomtypes':
        if len(line.strip().split()) == 0:
            fatomtypesMD.write(line)
            flag = 'NULL'
	else:
            fatomtypesMD.write(line)

    #Process the [ moleculestype ] section.
    elif line.startswith('[ moleculetype ]'):
        flag = 'moleculetype'
        ftopLMD.write(line)
    elif flag == 'moleculetype':
        if len(line.strip().split()) == 0:
            ftopLMD.write(line)
            flag = 'COPY'
	elif line.startswith(';'):
	    ftopLMD.write(line)
        else:
            ftopLMD.write("%s\t3\n" % nameL)

    #Copy the [ atoms ], [ bonds ], [ angles ], [ dihedrals ], and [ pairs ] sections. And write the restrait section.
    elif flag == 'COPY':
        if line.startswith('[ system ]'):
            ftopLMD.write("; Note the position restraint file can be inculded directly here, instead of the in the complex topology file.\n")
	    ftopLMD.write("; Include Position restraint file.\n")
	    ftopLMD.write("#ifdef POSRES\n")
	    ftopLMD.write("#include \"posre_LIG.itp\"\n")
	    ftopLMD.write("#endif\n")
	    break
        else:
	    ftopLMD.write(line)

ftopL0.close()
ftopLMD.close()
fatomtypesMD.close()

##########################################################################################
#Combine the toplogy of file of the ligand and the protein.
ftop = open(top,'r')
fatomtypesMD = open('LIG_atomtypes_MD.top','r')
ftopMD = open(topMD,'w')  #For classical MD

ftopMD.write("; Topology for classical MD.\n")

flag = 'NULL'
for line in ftop:
    #Incorporate the [ atomtypes ] section of the ligand after the force field section.
    if line.startswith('; Include forcefield parameters'):
        flag = 'include'
        ftopMD.write(line)
    elif flag == 'include': 
        if len(line.strip().split()) == 0: 
            ftopMD.write(line)
            for lineTMD in fatomtypesMD:
                ftopMD.write(lineTMD)
            flag = 'NULL'
        else:
            ftopMD.write(line)

    #Incorporate the ligand topology file
    elif line.startswith('; Include chain topologies'):
        flag = 'chains'
        #Include Ligand topology and restraints for ligand in Classical MD. This is just for the toplogy of proteins with multiple chains.
	#For proteins with multiple chains, it is more conveniet to put ligand topology before the protein chains.
        ftopMD.write("; Include ligand topology\n")
        ftopMD.write("#include \"%s_MD.top\"\n\n" % nameL)
	ftopMD.write(line)
    elif flag == 'chains':
        if len(line.strip().split()) == 0:
	    ftopMD.write(line)
            flag = 'COPY'
        else:
            ftopMD.write(line)

    #Copy the water and system section
    elif flag == 'COPY':
	if line.startswith('[ molecules ]'):
	    flag = 'molecules'
	    ftopMD.write(line)
	else:
	    ftopMD.write(line)

    #Process the [ molecules ] section    
    elif flag == 'molecules':
	if line.startswith(';'):
	    ftopMD.write(line)
	    ftopMD.write(nameL+"\t"+"1"+"\n")
	else:
            ftopMD.write(line)

ftop.close()
fatomtypesMD.close()
ftopMD.close()
##################################################################################
