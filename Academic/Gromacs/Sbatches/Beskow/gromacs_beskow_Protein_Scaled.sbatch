#!/bin/bash -l

## Name of the script
#SBATCH -J MD

## Names the error and output files according to the jobid
#SBATCH --error=EMNVTNPTMD_%J.err
#SBATCH --output=EMNVTNPTMD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=23:59:50                #Requested walltime 
#SBATCH --ntasks 128
#SBATCH --nodes=4 
#SBATCH --exclusive          #Number of compute nodes to allocate
#SBATCH --account=m.2017-1-425       #Tu, Account string for the project that you wish to account the job to
###SBATCH --account=2017-12-49     #Hans

#load the gromacs module
export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 128 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2016.5-plumed_2.3.5

###EM###
mkdir EM
aprun -n 1 gmx_mpi grompp -f em_protein.mdp -c complex_wi.gro -p complex.top -o EM/em.tpr
cd EM 
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm em
cd ..

###NVT###
mkdir NVT
aprun -n 1 gmx_mpi grompp -f nvt_protein.mdp -c EM/em.gro -p complex.top -n index.ndx -o NVT/nvt.tpr 
cd NVT
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm nvt
cd ..

###NPT###
mkdir NPT
aprun -n 1 gmx_mpi grompp -f npt_protein.mdp -c NVT/nvt.gro -t NVT/nvt.cpt -p complex.top -n index.ndx -o NPT/npt.tpr 
cd NPT
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm npt
cd ..

##Preprocess for scaling:
aprun -n 1 gmx_mpi grompp -f scaled_a7.mdp -c NPT/npt.gro -p complex.top -n index.ndx -pp complex_pp.top
/afs/pdc.kth.se/home/g/guanglin/nobackup/software/anaconda2/bin/python prepare_potential_scaled.py

###Scaled###
mkdir Scaled
aprun -n 1 gmx_mpi grompp -f scaled_a7.mdp -c NPT/npt.gro -t NPT/npt.cpt -p complex_scaled.top -n index.ndx -o Scaled/scaled.tpr
cp index.ndx gromacs_beskow_Protein_Scaled_Continue.sbatch Scaled
cd Scaled
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm scaled -plumed plumed_scaled.dat
cd ..

