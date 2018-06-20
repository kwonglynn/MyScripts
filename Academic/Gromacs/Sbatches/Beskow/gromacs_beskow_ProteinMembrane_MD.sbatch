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
module add gromacs/2016.3-plumed_2.3

###EM###
mkdir EM
aprun -n 1 gmx_mpi grompp -f em_membrane.mdp -c complex_wi.gro -p topol.top -o EM/em.tpr
cd EM 
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm em
cd ..

###NVT###
mkdir NVT
aprun -n 1 gmx_mpi grompp -f nvt_membrane.mdp -c EM/em.gro -p topol.top -n index.ndx -o NVT/nvt.tpr 
cd NVT
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm nvt
cd ..

###NPT###
mkdir NPT
aprun -n 1 gmx_mpi grompp -f npt_membrane.mdp -c NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT/npt.tpr 
cd NPT
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm npt
cd ..

###MD###
mkdir MD
aprun -n 1 gmx_mpi grompp -f md_membrane.mdp -c NPT/npt.gro -t NPT/npt.cpt -p topol.top -n index.ndx -o MD/md.tpr
cd MD
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm md
cd ..
