#1. Process the capping NME residue.
sed -e "1,\$s/CA  NMA/CH3 NMA/g" -e "1,\$s/NMA/NME/g" 2kkw-Helix-cap.pdb > 2kkw-Helix-cap-NME.pdb 

#2, Process the protein
gmx_mpi pdb2gmx -f 2kkw-Helix-cap-NME.pdb -o protein_processed.gro -ignh -ter

#3. Generate the box
# Make the box bigger because the peptide will unfold.
gmx_mpi editconf -f protein_processed.gro -o protein_box.gro -bt dodecahedron -d 2.0 

#4. Add water
gmx_mpi solvate -cp protein_box.gro -cs spc216.gro -p topol.top -o protein_w.gro

#5. Add ions
gmx_mpi grompp -f em_protein.mdp -c protein_w.gro -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o protein_wi.gro -p topol.top -neutral -conc 0.15

#7. Generate the index file for the system
echo -e "a OW\nq" | gmx_mpi make_ndx -f protein_wi.gro -o index.ndx

#10. Copy the input files for MD and Metadynamics simulations.
mkdir -p copy
cp -r charmm36-jul2017.ff protein_wi.gro topol.top posre*.itp index.ndx *mdp copy
cp -r copy ../Unbiased/
