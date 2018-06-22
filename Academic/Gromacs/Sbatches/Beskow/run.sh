for i in {1..20}; do cp -r copy/ Run$i; done
for i in {1..20}; do cd Run$i; sbatch gromacs_beskow_Protein_Scaled.sbatch; cd ..; done
for i in {1..20}; do cd Run$i; sbatch gromacs_beskow_Protein_Scaled_Continue.sbatch; cd ..; done
