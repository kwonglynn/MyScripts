for i in {1..20}; do cp -r copy/ Run$i; done
for i in {1..20}; do cd Run$i; sbatch gromacs_beskow_METAD.sbatch; cd ..; done
for i in {1..20}; do cd Run$i; sbatch gromacs_beskow_METAD_Restart.sbatch; cd ..; done

