for (( i=1;i<=20;i++))
do
    cd Node$i
    sbatch gromacs_kebnek_CPU_MD_${i}.sbatch
    cd ..
done 
