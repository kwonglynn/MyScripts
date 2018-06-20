for (( i=1;i<=20;i++))
do
    mkdir Node$i
    n=$(( i*14 ))
    sed -e "1,\$s/XXX/$i/g" -e "1,\$s/YYY/$n/g" gromacs_kebnek_CPU_MD.sbatch > Node${i}/gromacs_kebnek_CPU_MD_$i.sbatch
    cp md1.tpr Node$i
done 
