for (( i=4;i<=20;i++))
do
    mkdir Node$i
    n=$(( i*32 ))
    sed -e "1,\$s/XXX/$i/g" -e "1,\$s/YYY/$n/g" gromacs_beskow_MD1.sbatch > Node${i}/gromacs_beskow_MD1_$i.sbatch
    cp md1.tpr Node$i
done 
