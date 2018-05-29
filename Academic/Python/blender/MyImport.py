import imp
#Import the two libraries: Atoms.py and Molecule.py - change the path to where you keep them
atoms = imp.load_source('Atoms', './Atoms.py')
geometry = imp.load_source('Molecule', './Molecule.py')
import math

#path to the folder and file name for your structure
folder='./'
f='coordinates'

#my input filname is 'coordinates.xyz'
#this files contains the number of atoms on the first line;
#the 2nd line is blank
#then follow lines of thys type: symbol   x   y   z
fileName=folder+f+'.xyz'

#create a molecule object
molecule=geometry.Geometry()

#read the xyz file; Note: you can also read a vasp outcar file; the function is read_VASP(fileName, Nxyz=[na, nb, nc]
#where na, nb, nc define how many units cells you want to repeat in each direction
molecule.read_xyz(fileName)

#create molecule in blender; you can play with the scaling factor and bond size to change the dimensions of atoms and bonds;
#the keyword "threshold" defines the largest distance between two atoms that is considered a bond; it's in the units you use 
#in your .xyz file 
molecule.blender_draw_smart(scaling_factor=0.04,bond_size=0.2,diffCol=[],threshold=1.7,scale=1.0)
	
#Uncomment below if you want to write a ps file that can be imported in inkscape:
###Translate the molecule to the center of the page (assumes that your system is centered at 0)
##molecule.translate_along(12.5, 12.5, 12.5)
###Write ps file
##Name=folder+f+'.ps'
##molecule.write_ps_smart(psName=Name, step=25)


