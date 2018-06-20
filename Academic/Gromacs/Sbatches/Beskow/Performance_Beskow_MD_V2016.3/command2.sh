for (( i=4;i<=20;i++))
do
    cd Node$i
    sbatch gromacs_beskow_MD1_${i}.sbatch
    cd ..
done 
