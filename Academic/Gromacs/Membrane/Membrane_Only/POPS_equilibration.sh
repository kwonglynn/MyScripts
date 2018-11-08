#1, Generate the coordinates of the POPS membrane with CHRAMMGUI, the x, y length is 8 nm

#2, Convert the PDB (step5_charmm2gmx.pdb) to gromacs
grep POPS step5_charmm2gmx.pdb > POPS_GUI.pdb
gmx_mpi editconf -f POPS_GUI.pdb -o POPS_GUI.gro

#3. Measure the minmax of POPS_GUI.gro, get the box size for POPS, and measure the center of POPS_GUI.gro to get the center for the box

#4. Generate the box
gmx_mpi editconf -f POPS_GUI.gro -o POPS_box.gro -box 8.0 8.0 8.0 -center 4.0 4.0 4.0

#5. Add water
gmx_mpi solvate -cp POPS_box.gro -cs spc216.gro -p topol.top -o POPS_w.gro

#6. View the structure in VMD, change the size of the box if necessary.

#7. Rename or remove vdwradii.dat.
mv vdwradii.dat KGL_vdwradii.dat

#8. Add ions
gmx_mpi grompp -f em_membrane.mdp -c POPS_w.gro -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o POPS_wi.gro -p topol.top -neutral -conc 0.15

#9. Extract the head part of POPS in pymol and save it as a PDB file. Then extract the atom names of the head atoms.

#10. Make an index file for the head atoms
gmx_mpi make_ndx -f POPS-1_processed.gro -o index-POPS-1.ndx
a N P C1 C2 C3 C11 O11 C12 O12 C13 O13 O13A O13B O14 O21 O31

#11. Make a restraint file for the head atoms
gmx_mpi genrestr -f POPS-1_processed.gro -o posre_Head.itp -n index-POPS-1.ndx -fc 1000 1000 1000

#12. Incorporate the retraint file for the head atoms in the topology file of POPS, namely POPS.top.

#13. Make a index file for the membrane-water system.
gmx_mpi make_ndx -f POPS_wi.gro -o index.ndx

#14. Create two mdp files for density equilibration.
#In npt1_membrane.mdp, all the membrane atoms are restrained, for 1 ns.
#In npt2_membrane.mdp, only the head atoms are restrainted, for 1 ns.

#15. In production run, no restraint is applied, for 100 ns.
