#1. Generate the topology file for the protein.
# gmx_mpi pdb2gmx -f 2N0A-KGL.pdb -o 2N0A_processed.gro -ignh -ter

#Part I: Apo
#1. Generate the box
gmx_mpi editconf -f 2N0A_processed.gro -o protein_box.gro -bt dodecahedron -d 1.0
#2. Add water
gmx_mpi solvate -cp protein_box.gro -cs spc216.gro -p topol.top -o protein_w.gro
#3. Add ions
gmx_mpi grompp -f em_protein.mdp -c protein_w.gro -p topol.top -o ions.tpr
echo 13 | gmx_mpi genion -s ions.tpr -o protein_wi.gro -p topol.top -neutral -conc 0.15
#4. Generate the restraint file for the protein backbone
echo 4 | gmx_mpi genrestr -f protein_wi.gro -o posre_backbone_raw.itp -fc 1000 1000 1000
#Delete the lines of posre_backone_raw.itp after number 890, the remaining is for chain A, and it should be noted that every chain is exactly the same.
mv posre_backbone_raw.itp posre_Protein_chain_A_backbone.itp
cp posre_Protein_chain_A_backbone.itp posre_Protein_chain_B_backbone.itp
...
#Add the following lines in topol_Protein_Chain_A.itp, and for other chains in the same way.

; Include Position restraint file
#ifdef POSRES_backbone
#include "posre_Protein_chain_A_backbone.itp"
#endif

#Use the script to do it in batch.
for i in B C D E F G H I J; do cp posre_Protein_chain_A_backbone.itp posre_Protein_chain_${i}_backbone.itp; done
for i in A B C D E F G H I J; do sed -e "1,\$s/chain_A/chain_${i}/g" posre_backbone.txt > posre_backbone_${i}.txt; done
for i in A B C D E F G H I J; do cat posre_backbone_${i}.txt >> topol_Protein_chain_${i}.itp; done

#Part II: Holo
#2. Use the Process_MD_Multiple_V1.py to process the topology and coordinate files.
# python Process_MD_Multiple_V1.py
#3. Generate the box
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
