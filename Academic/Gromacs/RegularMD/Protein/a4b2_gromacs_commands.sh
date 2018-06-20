#1. PDB to GMX
gmx_mpi pdb2gmx -f protein.pdb -o 5KXI_processed.gro -ignh -ter
#1. Generate the box
gmx_mpi editconf -f 5KXI_2FA.gro -o 5KXI_2FA_box.gro -bt dodecahedron -d 1.0
#2. Add water
gmx_mpi solvate -cp 5KXI_2FA_box.gro -cs spc216.gro -p topol.top -o 5KXI_2FA_w.gro
#3. Add ions
gmx_mpi grompp -f em_protein.mdp -c 5KXI_2FA_w.gro -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o 5KXI_2FA_wi.gro -p topol.top -neutral -conc 0.15
#4. Generate the restraint file for the ligand
gmx_mpi genrestr -f 5KXI_2FA_wi.gro -o posre_LIG.itp -fc 1000 1000 1000
#5. Generate the index file for the system
gmx_mpi make_ndx -f 5KXI_2FA_wi.gro -o index.ndx
