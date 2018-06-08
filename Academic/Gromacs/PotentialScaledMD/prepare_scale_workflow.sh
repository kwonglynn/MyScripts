##1. Generate an mdp file potential scaled MD.
#Define the restraints for all the protein residues, except those within 5 of angstroms of the ligand.
#Use NVT for scaled MD.

##2. Define the atoms to be restrained, exluding the atoms within 5 angstrom of the ligand.
#The flexible protein atoms are located on two subsequnt chains.
#Delete the lines for the flexible atoms in the relevant posre_ files.
#Do it for the larget ligand KI-n89, once and for all.

##3. Generate a topology file containing all the force field paremeters with grompp -pp
gmx_mpi grompp -f scaled_a7.mdp -c NPT/npt.gro -p complex.top -n index.ndx -pp processed.top
