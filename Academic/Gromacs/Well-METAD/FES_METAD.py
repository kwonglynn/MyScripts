###############Input####################################
#Command for plumed
#plumed sum_hills --hills HILLS --bin 98,99
fin='fes.dat'
fout='fes1.png'
auto=0      #0 for use input, 1 for automatic generation.
upper=1     #Always be 1
lower=-17   #1 smaller than the automatic one
STEP=1
bin=99
##########Input Stop here###############################

import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.tri as tri
import numpy as np
import math

data = np.genfromtxt(fin)
allX=data[:,0]*10
allY=data[:,1]
allZ=data[:,2]/4.184 		#kJ/mol->kcal/mol

zmax=allZ.max()
allZ0=allZ-zmax			#Let the upper limit to be 0

if auto==1:
    lower=int(allZ0.min())
    upper=int(allZ0.max())

X=allX[0:bin]
Y=allY[range(1,len(allY),bin)]
Z=allZ0.reshape(bin,bin)

levels=range(lower,upper,STEP)
print "levels are:\n",levels

CS=plt.contourf(X,Y,Z,levels,cmap='jet',origin='lower')

cbar=plt.colorbar(CS)

xmin=X.min()
xmax=X.max()
ymin=Y.min()
ymax=Y.max()

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.xlabel('Distance' +'('+r'$\AA$'+')',fontsize=20)
plt.ylabel('Distance' +'('+r'$\AA$'+')',fontsize=20)
cbar.ax.tick_params(labelsize=18)
cbar.set_label(r'$\Delta$'+'G'+' (kcal/mol)',fontsize=18)
plt.tick_params(labelsize=18)
plt.savefig(fout, dpi=600)

