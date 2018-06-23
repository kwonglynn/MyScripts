for i in {0..29}: do
    if [ -f complex_$i.gro ]; then
        mkdir Umbrella$i
        cp complex_${i}.gro common/* Umbrella$i
        cd Umbrella$i
        cp complex${i}.gro complex.gro
        sh gromacs_common.sh
        echo '[ ref ]' >> index.ndx
        echo '4734 4750 4760 4774 4815 4837 4851 4873 5021 5035 5042' >> index.ndx
        cd ..
    fi
done

for i in {0..29}; do
    if [ -f complex_$i.gro ]; then
        cd Umbrella${i}
        sbatch gromacs_beskow_Protein_Umbrella.sbatch
        cd ..
    fi
done


for i in {0..29}; do
    if [ -f complex_$i.gro ]; then
    cd Umbrella${i}/Umbrella
    cp umbrella.tpr umbrella${i}.tpr
    cp umbrella_pullf.xvg umbralla_pullf${i}.xvg
    cp umbrella${i}.tpr umbralla_pullf${i}.xvg ../../Analysis
    cd ../..
    fi
done

cd Analysis
ls umbrella*tpr > tpr-files.dat
ls umbrella_pullf*xvg > pullf-filex.dat
gmx_seq wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
xmgrace -nxy histogram.xvg