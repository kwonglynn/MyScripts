"""
Author:
    Guanglin Kuang <guanglin@kth.se>
    JUN 12, 2018

Usage:
    Preprocess_MD_Multiple_V2.py [(-i <insert> --num_before_LIG <num_before_LIG>)] [options]
    Preprocess_MD_Multiple_V2.py -h | --help

Options:
    -i, --insert <insert>                   The insert position of the ligand, see notes below [default: 1].
    --num_before_LIG <num_before_LIG>       The number of protein atoms before the ligand [default: 0].
    --coord_prot <coord_prot>               The coordinates of the protein [default: protein_processed.gro].
    --coord_LIG <coord_LIG>                 The coordinates of the ligand [default: LIG_MD.gro].
    --coord_com <coord_com>                 The coordinates of the complex [default: complex.gro].
    --top_prot <top_prot>                   The topology file of the protein [default: topol.top].
    --top_LIG0 <top_LIG0>                   The raw topology file of the ligand [default: LIG_GMX.top].
    --top_LIG <top_LIG>                     The processed topology file of the ligand [default: LIG_MD.top].
    --top_com <top_com>                     The topology file of the complex [default: complex.top].
    -h, --help                              Show this screen.

Notes:
    If put the ligand before all the chains, insert = 1, namely, the ligand is the first chain.
    If put the ligand after the first chain, insert = 2, namely, the ligand in the second chain.
    Etc...

Common:
    python Preprocess_MD_Multiple_V2.py -i 1 --num_before_LIG 0 --coord_prot protein_processed.gro --coord_LIG LIG_MD.gro --top_prot topol.top --top_LIG LIG_GMX.top
    
"""

from docopt import docopt

#Combine the coordinate files of protein and ligand.
def combine_coord(coord_prot, coord_LIG, coord_com, num_before_LIG):
    fi_coordP = open(coord_prot, 'r')
    fi_coordL = open(coord_LIG, 'r')
    fo_coordC = open(coord_com, 'w')
    
    xyz_prot = fi_coordP.readlines()
    xyz_LIG = fi_coordL.readlines()
    
    NP = len(xyz_prot) - 3
    NL = len(xyz_LIG) - 3
    NC = NP + NL
    
    fo_coordC.write("Gromacs coordinates for protein and ligand complex\n")
    fo_coordC.write(str(NC) + "\n")
    
    fo_coordC.writelines(xyz_prot[2:-1][:num_before_LIG])
    fo_coordC.writelines(xyz_LIG[2:-1])
    fo_coordC.writelines(xyz_prot[2:-1][num_before_LIG:])
   
    fo_coordC.write(xyz_prot[-1])
    
    fi_coordP.close()
    fi_coordL.close()
    fo_coordC.close()

#Process the ligand topology file.
#The blank lines of the topology files are very important! Don't change the default settings!
def process_LIG(top_LIG0, top_LIG):
    fi_topL0 = open(top_LIG0, 'r')  #The "raw" topology file of the ligand produced by antechamber and acpype.
    fo_topL = open(top_LIG, 'w')
    fo_atomtypes = open('LIG_atomtypes_MD.top','w')
    
    flag = 'NULL'
    for line in fi_topL0:
        # Judge the section type.
        if line.startswith('[ atomtypes ]'):
            flag = 'atomtypes'
            # Note that the [ atomtype ] section is written to the complex topology, not the ligand topology.
            fo_atomtypes.write(line)    
        elif line.startswith('[ moleculetype ]'): 
            # Copy from [ moleculetypes ] until [ system ] sections.
            flag = 'COPY'
            fo_topL.write(line)            
        elif line.startswith('[ system ]'):
            flag = 'END'

        # Process corresponding sections:
        #Process the [ atomtypes ] section.        
        elif flag == 'atomtypes':
            fo_atomtypes.write(line)
        
        #Process the [ moleculestype ] section.
        elif flag == 'COPY':
            fo_topL.write(line)
    
        #Copy the [ atoms ], [ bonds ], [ angles ], [ dihedrals ], and [ pairs ] sections. And write the restrait section.
        elif flag == 'END':
            restraint = '''
; Note the position restraint file can be inculded directly here, instead of the in the complex topology file.
; Include Position restraint file.
#ifdef POSRES
#include "posre_LIG.itp"
#endif\n
'''
            fo_topL.write(restraint)
            break
    
    fi_topL0.close()
    fo_atomtypes.close()
    fo_topL.close()
    
#Combine the toplogy of file of the ligand and the protein.
def combine_top(top_Prot, top_com, insert):
    fi_topP = open(top_Prot, 'r')
    fi_atomtypes = open('LIG_atomtypes_MD.top', 'r')
    fo_topC = open(top_com, 'w')  #For unbiased MD
    
    fo_topC.write("; Topology for unbiased MD.\n")
    
    lines = fi_topP.readlines()
    #for line in ftop:
    i = 0
    while i < len(lines):
        if lines[i].startswith('; Include forcefield parameters'):
            fo_topC.writelines(lines[i:i+3])
            for line_at in fi_atomtypes:
                fo_topC.write(line_at)
            i = i + 3
            
        elif lines[i].startswith('; Include chain topologies'):
            fo_topC.writelines(lines[i:i+insert])
            fo_topC.write("\n; Include ligand topology\n")
            fo_topC.write("#include \"LIG_MD.top\"\n\n")
            i = i + insert
        
        elif lines[i].startswith('[ molecules ]'):
            fo_topC.writelines(lines[i:i+1+insert]) #There is a comment line before the protein.           
            fo_topC.write('LIG' + "\t" + "1" + "\n")
            i = i + 1 + insert
        
        else:
            fo_topC.write(lines[i])
            i = i + 1
            
    fi_topP.close()
    fi_atomtypes.close()
    fo_topC.close()
##################################################################################
if __name__ == '__main__':
    opts = docopt(__doc__)
    
    insert = int(opts["--insert"])
    num_before_LIG = int(opts["--num_before_LIG"])
    coord_prot = opts["--coord_prot"]
    coord_LIG = opts["--coord_LIG"]
    coord_com = opts["--coord_com"]
    top_prot = opts["--top_prot"]
    top_LIG0 = opts["--top_LIG0"]
    top_LIG = opts["--top_LIG"]
    top_com = opts["--top_com"]

    if coord_prot:    
        combine_coord(coord_prot, coord_LIG, coord_com, num_before_LIG)
    
    if top_LIG0:
        process_LIG(top_LIG0, top_LIG)
    
    if top_prot:
        combine_top(top_prot, top_com, insert)