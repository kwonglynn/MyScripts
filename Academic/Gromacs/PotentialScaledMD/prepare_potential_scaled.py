# -*- coding: utf-8 -*-

"""
Author:
    Guanglin Kuang <guanglin@kth.se>
    JUN 10, 2018

Usage:
    process_potential_scaled.py [options]

Options:
    -i, --input <input>         The toplogy file of the system processed with -pp [default: complex_pp.top].
    -l, --lambda <lambda>       The scaling factor lambda [default: 0.4].
    -o, --output <file>         The output file. [default: complex_scaled.top].

Note:
    Command to generate processed.top:
    gmx_mpi grompp -f scaled.mdp -c NPT/npt.gro -p complex.top -n index.ndx -pp Scaled/complex_pp.top
    Potential scaling is only used in the production run with NVT ensemble, therefore the densitiy has to be equilibrated first.

"""
import numpy as np
from docopt import docopt

opts = docopt(__doc__)

fi = open(opts["--input"], 'r')
fo = open(opts["--output"], 'w')
scale = float(opts["--lambda"])     # Scaling factor.

flag = 'NULL'
for line in fi:
    # Lines starting with '[' are directives. Read the line, allocate a flag for further processing and write the line.
    # Process the lines from top to bottom, one by one. Always remember what kind of line is being processed.
    if line.startswith('['):
        fo.write(line)
        if line.startswith('[ atomtypes ]'):
            flag = 'atomtypes'
        elif line.startswith('[ bondtypes ]'):
            flag = 'bondtypes'
        elif line.startswith('[ angletypes ]'):
            flag = 'angletypes'            
        elif line.startswith('[ dihedraltypes ]'):
            flag = 'dihedraltypes'
        elif line.startswith('[ atoms ]'):
            flag = 'atoms'
        elif line.startswith('[ bonds ]'):
            flag = 'bonds'           
        elif line.startswith('[ angles ]'):
            flag = 'angles'
        elif line.startswith('[ dihedrals ]'):
            flag = 'dihedrals'
        else:
            #Don't forget other directives that don't need to be scaled.
            flag = 'ELSE'

    # Write comment lines or empty lines, ignore the flag!    
    elif line.startswith(';') or line.startswith('*') or len(line.strip()) == 0:
        fo.write(line)

    # Write Other directies that don't need to be scaled.
    elif flag == 'ELSE':
        fo.write(line)
    
    # For atomstypes, scale epsilon which is in column 6 (starting from 0). 
    # Note the return '\n' at the end of line is removed after splitting.
    elif flag == 'atomtypes':
        items = line.split()
        epsilon = float(items[6]) * scale
        newitems = items[:6] + ["%.5e" % epsilon] 
        newline = '\t'.join(newitems) + '\n'
        fo.write(newline)

    # For bondtypes, scale kb which is in cloumn 4.      
    elif flag == 'bondtypes':
        items = line.split()
        kb = float(items[4]) * scale
        newitems = items[:4] + ["%.1f" % kb] + items[5:]
        newline = '\t'.join(newitems) + '\n'
        fo.write(newline)

    # For angletyples, scale cth which is in column 5.
    elif flag == 'angletypes':
        items = line.split()
        cth = float(items[5]) * scale
        newitems = items[:5] + ["%.3f" % cth] + items[6:]
        newline = '\t'.join(newitems) + '\n'
        fo.write(newline)

    # For dihedraltypes, scale kd which in column 6.
    # This is for the dihedraltypes of the protein. The dihedrals for ligands are written explicitly in the [ dihedrals] section.
    # For protein dihedrals, all the parameters that need to be scaled are located at column 6, even if with different function types..
    elif flag == 'dihedraltypes':
        items = line.split()
        kd = float(items[6]) * scale
        newitems = items[:6] + ["%.5f" % kd] + items[7:]
        newline = '\t'.join(newitems) + '\n'
        fo.write(newline)

    # For atoms, scale charge which is in cloumn 6 with the squred root of lambda.        
    elif flag == 'atoms':
        items = line.split()
        charge = float(items[6]) * np.sqrt(scale)
        newitems = items[:6] + ["%.6f" % charge] + items[7:]
        newline = '\t'.join(newitems) + '\n' 
        fo.write(newline)                   

    # The bonds for protein and ligand are processed separately.             
    elif flag == 'bonds':
        # For protein, no explicit parameters for bonds are written.
        items = line.split()
        if len(items) == 3:
            fo.write(line)
        elif len(items) > 3:
            # For ligands, scale k which is in column 4.
            k = float(items[4]) * scale
            newitems = items[:4] + ["%.4e" % k] + items[5:]
            newline = '\t'.join(newitems) + '\n'
            fo.write(newline)

    # The angles for protein and ligand are processed separately.             
    elif flag == 'angles':
        # For protein, no explicit parameters for angles are written.
        items = line.split()
        if len(items) == 4:
            fo.write(line)
        elif len(items) > 4:
            # For ligands, scale cth which is in column 5.
            cth = float(items[5]) * scale
            newitems = items[:5] + ["%.4e" % cth] + items[6:]
            newline = '\t'.join(newitems) + '\n'
            fo.write(newline)

    # The dihedrals for protein and ligand are processed separately.                 
    elif flag == 'dihedrals':
        # For protein (except ILDN), no explicit parameters for dihedrals are written.
        items = line.split()
        if len(items) == 5:
            fo.write(line)
        elif len(items) > 5:
        # For the ligand, the dihedrals include proper and improper ones.
            func = int(items[4])
            # For proper dihedrals, scale C0 to C5 (column 5 to 10)
            if func == 3:
                c0 = float(items[5]) * scale
                c1 = float(items[6]) * scale
                c2 = float(items[7]) * scale
                c3 = float(items[8]) * scale
                c4 = float(items[9]) * scale
                c5 = float(items[10]) * scale
                newitems = items[:5] + ["%.5f" % c0] + ["%.5f" % c1] + ["%.5f" % c2] + ["%.5f" % c3] + ["%.5f" % c4] + ["%.5f" % c5] + items[11:]
                newline = '\t'.join(newitems) + '\n'
                fo.write(newline)
            # For improper dihedrals, scale kd which is in column 6
            # func == 9 is the ILDN terms for protein.
            elif func == 1 or func == 9:
                kd =float(items[6]) * scale
                newitems = items[:6] + ["%.5f" % kd] + items[7:]
                newline = '\t'.join(newitems) + '\n'              
                fo.write(newline)
                
fi.close()
fo.close()    