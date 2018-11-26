#!/bin/bash
#SBATCH -A SNIC2017-12-49
#SBATCH -J Amber-GPU
#SBATCH -n 2
#SBATCH --gres=gpu:k80:1
#SBATCH --time=7-00:00:00

# It is always best to do a ml purge before loading modules in a submit file
ml purge > /dev/null 2>&1
module load GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
module load Amber/18-AmberTools-18-patchlevel-10-8

#(1) Minimization with strong restraints.
export CUDA_VISIBLE_DEVICES=0
cd Folder1
pmemd.cuda -O -i 01_Min.in -o 01_Min.out -p complex_wi.prmtop  -c complex_wi.inpcrd -r 01_Min.rst -ref complex_wi.inpcrd

#(2) Minimization with no restraints.
pmemd.cuda -O -i 02_Min.in -o 02_Min.out -p complex_wi.prmtop  -c 01_Min.rst -r 02_Min.rst

#(3) Heating with weak restraints.
pmemd.cuda -O -i 03_Heat.in -o 03_Heat.out -p complex_wi.prmtop -c 02_Min.rst -r 03_Heat.rst -ref 02_Min.rst -x 03_Heat.nc

#(4) Density Equilibration with weak restraints.
pmemd.cuda -O -i 04_Density.in -o 04_Density.out -p complex_wi.prmtop -c 03_Heat.rst -r 04_Density.rst -ref 03_Heat.rst -x 04_Density.nc

#(5) Production with no restraints.
pmemd.cuda -O -i 05_Prod.in -o 05_Prod.out -p complex_wi.prmtop -c 04_Density.rst -r 05_Prod.rst -x 05_Prod.nc &

export CUDA_VISIBLE_DEVICES=1
cd Folder2
pmemd.cuda -O -i 01_Min.in -o 01_Min.out -p complex_wi.prmtop  -c complex_wi.inpcrd -r 01_Min.rst -ref complex_wi.inpcrd

#(2) Minimization with no restraints.
pmemd.cuda -O -i 02_Min.in -o 02_Min.out -p complex_wi.prmtop  -c 01_Min.rst -r 02_Min.rst

#(3) Heating with weak restraints.
pmemd.cuda -O -i 03_Heat.in -o 03_Heat.out -p complex_wi.prmtop -c 02_Min.rst -r 03_Heat.rst -ref 02_Min.rst -x 03_Heat.nc

#(4) Density Equilibration with weak restraints.
pmemd.cuda -O -i 04_Density.in -o 04_Density.out -p complex_wi.prmtop -c 03_Heat.rst -r 04_Density.rst -ref 03_Heat.rst -x 04_Density.nc

#(5) Production with no restraints.
pmemd.cuda -O -i 05_Prod.in -o 05_Prod.out -p complex_wi.prmtop -c 04_Density.rst -r 05_Prod.rst -x 05_Prod.nc &

wait
