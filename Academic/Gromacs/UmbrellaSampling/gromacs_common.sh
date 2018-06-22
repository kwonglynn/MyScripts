#1. Generate the water box
gmx_seq editconf -f complex.gro -o complex_box.gro -bt dodecahedron -d 1.0
#2. Add water
gmx_seq solvate -cp complex_box.gro -cs spc216.gro -p complex.top -o complex_w.gro
#3. Add ions
gmx_seq grompp -f em_protein.mdp -c complex_w.gro -p complex.top -o ions.tpr
echo 15 | gmx_seq genion -s ions.tpr -o complex_wi.gro -p complex.top -neutral -conc 0.15
#4. Generate the index file for the system
echo -e "2|5\na OW\nq" | gmx_seq make_ndx -f complex_wi.gro -o index.ndx
