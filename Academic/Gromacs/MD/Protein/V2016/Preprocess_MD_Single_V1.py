"""
Author:
    Guanglin Kuang <guanglin@kth.se>
    JUN 13, 2018

Usage:
    Preprocess_MD_Single_V1.py [options]
    Preprocess_MD_Single_V1.py -h | --help

Options:
    --coord_prot <coord_prot>               The coordinates of the protein [default: protein_processed.gro].
    --coord_LIG <coord_LIG>                 The coordinates of the ligand [default: LIG_MD.gro].
    --coord_com <coord_com>                 The coordinates of the complex [default: complex.gro].
    --top_prot <top_prot>                   The topology file of the protein [default: topol.top].
    --top_LIG0 <top_LIG0>                   The raw topology file of the ligand [default: LIG_GMX.top].
    --top_LIG <top_LIG>                     The processed topology file of the ligand [default: LIG_MD.top].
    --top_com <top_com>                     The topology file of the complex [default: complex.top].
    -h, --help                              Show this screen.

Common:
    python Preprocess_MD_Single_V1.py --before_protein --coord_prot protein_processed.gro --coord_LIG LIG_MD.gro --top_prot topol.top --top_LIG0 LIG_GMX.top
    
"""

from docopt import docopt

#Combine the coordinate files of protein and ligand.
def combine_coord(coord_prot, coord_LIG, coord_com):
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
    
    fo_coordC.writelines(xyz_prot[2:-1])        
    fo_coordC.writelines(xyz_LIG[2:-1])
        
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
def combine_top(top_Prot, top_com):
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
            
        elif lines[i].startswith('; Include water topology'):
            fo_topC.write("; Include ligand topology\n")
            fo_topC.write("#include \"LIG_MD.top\"\n\n")
            fo_topC.write(lines[i])
            i = i + 1                          
        
        elif lines[i].startswith('[ molecules ]'):
            fo_topC.writelines(lines[i:i+3]) # Write comment and the protein line first. 
            fo_topC.write('LIG' + "\t" + "1" + "\n")
            i = i + 3                
        
        else:
            fo_topC.write(lines[i])
            i = i + 1
            
    fi_topP.close()
    fi_atomtypes.close()
    fo_topC.close()
##################################################################################
if __name__ == '__main__':
    opts = docopt(__doc__)

    coord_prot = opts["--coord_prot"]
    coord_LIG = opts["--coord_LIG"]
    coord_com = opts["--coord_com"]
    top_prot = opts["--top_prot"]
    top_LIG0 = opts["--top_LIG0"]
    top_LIG = opts["--top_LIG"]
    top_com = opts["--top_com"]
    
    if coord_prot:    
        combine_coord(coord_prot, coord_LIG, coord_com)
    
    if top_LIG0:
        process_LIG(top_LIG0, top_LIG)
    
    if top_prot:
        combine_top(top_prot, top_com)