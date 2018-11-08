#1. Align the protein along the principal axis
gmx_mpi editconf -f protein_processed.gro -o protein_princ.gro -princ -c

#2, Generate the coordinates of the POPS membranes with CHRAMMGUI, the x, y length is 6 nm

#3, Convert the PDB (step5_charmm2gmx.pdb) to gromacs
grep POPS step5_charmm2gmx.pdb > POPS_GUI.pdb
gmx_mpi editconf -f POPS_GUI.pdb -o POPS_GUI.gro

#4, Extract one lipid PDB coordinate, generate its topology
grep 'POPS    1' POPS_GUI.pdb > POPS-1.pdb
gmx_mpi pdb2gmx -f POPS-1.pdb -o POPS-1_processed.gro -p POPS.top -i posre_POPS.itp
# Delete the unnecessary lines of POPS.top

#5, Move the protein to the top of the membrane, with the basic residues facing the membrane.
# Use VMD
soure move_protein.tcl

#6, Convert the moved proten to gromacs format.
gmx_mpi editconf -f protein_P2.pdb -o protein_P2.gro

#7, Combine the coordinates of the protein and the membrane, remove unnecessary lines and change the atom number on the second line.
cat protein_P2.gro POPS-noPBC-100ns-noWater.gro > complex.gro

#3. Generate the box
# Make the box bigger because the peptide will unfold.
gmx_mpi editconf -f complex.gro -o complex_box.gro -box 7.63 7.63 13 -center 3.81 3.81 4.50

#4. Add water
# Note the vdwradii.dat file is needed for adding water and ions, change the value of C from 0.15 to 0.375
gmx_mpi solvate -cp complex_box.gro -cs spc216.gro -p topol.top -o complex_w.gro

#5. Add ions
gmx_mpi grompp -f em_membrane.mdp -c complex_w.gro -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o complex_wi.gro -p topol.top -neutral -conc 0.15
# Note, the vdwradii.dat should be re-named, change it to Lipid-radii.dat

#7. Generate the index file for the system
echo -e "1|13\na OW\nq" | gmx_mpi make_ndx -f complex_wi.gro -o index.ndx

#10. Copy the input files for MD and Metadynamics simulations.
mkdir -p copy
cp -r charmm36-jul2017.ff complex_wi.gro topol.top posre*.itp index.ndx POPS.top plumed*.dat *mdp gromacs*.sbatch copy
cp -r copy ../Unbiased/
