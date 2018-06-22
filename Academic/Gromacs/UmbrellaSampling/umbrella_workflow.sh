for i in {1..20}: do
    mkdir Umbrella$i
    cp complex${i}.gro common/* Umbrella$i
    cd Umbrella$i
    cp complex${i}.gro complex.gro
    sh gromacs_common.sh
    cd ..
done

for i in {1..20}; do cd Umbrella${i}; sbatch gromacs_beskow_Protein_Umbrella.sbatch; done
for i in {1..20}; do 
    cd Umbrella${i}/Umbrella
    cp umbrella.tpr umbrella${i}.tpr
    cp pullf-umbrella.xvg pullf-umbralla${i}.xvg
    cp umbrella${i}.tpr pullf-umbralla${i}.xvg ../../Analysis
    cd ..
done

cd Analysis
ls umbrella*tpr > tpr-files.dat
ls pullf*xvg > pullf-filex.dat
gmx_seq wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
xmgrace -nxy histogram.xvg
