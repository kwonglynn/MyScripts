#!/usr/bin/env python
# encoding: utf-8

"""
Author:
    Guanglin Kuang <guanglin@kth.se>

Usage:
    FES_METAD.py <fes> [(-u <upper> -l <lower>)] --cv1 <cv1> --cv2 <cv2> [options]

Options:
    -o, --output <file>         Save the plot to a file [default: fes.png].
    -u, --upper <upper>         The upper boundary of the map.
    -l, --lower <lower>         The lower boundary of the map.
    -s, --step <step>           The step for the boundary [default: 1].
    -B, --BIN <BIN>             The bin number for the map [default: 99].
    --cv1 <cv1>                 CV type, can be distance, dihedral or contact [default: distance].
    --cv2 <cv2>                 CV type, can be distance, dihedral or contact [default: dihedral].
    
    -v, --verbose               Verbose mode.

Note:
    1, Use plumed to generate the free energy data:
       plumed sum_hills --hills HILLS --bin 98,99
       Distance BIN is 98, dihedral BIN is 99. Try and check fes.dat for the details.
    2, Use the script to generate a FES map.
    3, Adjust the upper and lower boundaries to make the map nicer.
"""

import matplotlib.pyplot as plt
import numpy as np
from docopt import docopt

##########Input Stop here###############################
if __name__ == '__main__':
    
    ####Treat the options####
    opts = docopt(__doc__)
    
    fin = opts["<fes>"]
    fout = opts["--output"]
    
    if opts["--upper"] and opts["--lower"]:
        upper = int(opts["--upper"])
        lower = int(opts["--lower"])
    
    step = int(opts["--step"])
    BIN = int(opts["--BIN"])
    #########Options#########
    
    data = np.genfromtxt(fin)
    
    # The first collective variable.
    if opts["--cv1"] == 'distance':
        allX = data[:, 0] * 10       # nM to Angstrom
        unit1 = r'$\AA$'
    elif opts["--cv1"] == 'dihedral':
        allX = data[:, 0]
        unit1 = 'radian'
    elif opts["--cv1"] == 'contact':
        allX = data[:, 0]
        unit1 = ''
    
    # The second collective variable.
    if opts["--cv2"] == 'distance':
        allY = data[:, 1] * 10       # nM to Angstrom
        unit2 = r'$\AA$'
    elif opts["--cv2"] == 'dihedral':
        allY = data[:, 1]
        unit2 = 'radian'
    elif opts["--cv2"] == 'contact':
        allY = data[:, 1]
        unit2 = ''        
        
    allZ = data[:, 2] / 4.184 	  # kJ/mol->kcal/mol
    
    zmax = allZ.max()
    allZ0 = allZ - zmax			#Let the upper limit to be 0

    # Either provide both upper and lower boundaries or leave it empty.
    if not opts["--upper"]:
        lower0 = int(allZ0.min())
        upper0 = int(allZ0.max())
    
        lower = lower0 - 1
        upper = upper0 + 1
    
    X = allX[0 : BIN]
    Y = allY[range(1, len(allY), BIN)]
    Z = allZ0.reshape(BIN, BIN)
    
    levels = range(lower, upper, step)
    print ("Levels are:\n", levels)
    
    CS = plt.contourf(X, Y, Z, levels, cmap='jet', origin='lower')
    
    cbar = plt.colorbar(CS)
    
    xmin = X.min()
    xmax = X.max()
    ymin = Y.min()
    ymax = Y.max()
    
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    
    # Customize the x, y label if necessary to make it clearer.
    plt.xlabel(opts["--cv1"].capitalize() + '(' + unit1 + ')', fontsize=12)
    plt.ylabel(opts["--cv2"].capitalize() + '(' + unit2 + ')', fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(r'$\Delta$' + 'G' + ' (kcal/mol)', fontsize=10)
    plt.tick_params(labelsize=10)
    plt.savefig(fout, dpi=600)
    
