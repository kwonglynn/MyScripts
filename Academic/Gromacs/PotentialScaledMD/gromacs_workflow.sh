gmx_mpi editconf -f complex.gro -o complex_box.gro -bt dodecahedron -d 1.0
#4. Add water
gmx_mpi solvate -cp complex_box.gro -cs spc216.gro -p complex.top -o complex_w.gro
#5. Add ions
gmx_mpi grompp -f em_protein.mdp -c complex_w.gro -p complex.top -o ions.tpr
echo 15 | gmx_mpi genion -s ions.tpr -o complex_wi.gro -p complex.top -neutral -conc 0.15
#6. Generate the restraint file for the ligand
echo 2 | gmx_mpi genrestr -f complex_wi.gro -o posre_LIG.itp -fc 1000 1000 1000
#7. Generate the index file for the system
echo -e "2|5\na OW\nq" | gmx_mpi make_ndx -f complex_wi.gro -o index.ndx
mkdir -p copy
cp LIG_MD.top complex.top complex_wi.gro index.ndx *mdp posre_* topol* copy/
