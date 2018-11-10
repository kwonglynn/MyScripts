#!/bin/bash -l

## Name of the script
#SBATCH -J REST2

## Names the error and output files according to the jobid
#SBATCH --error=REST2_%J.err
#SBATCH --output=REST2_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=23:59:50                #Requested walltime 
#SBATCH --ntasks 320
#SBATCH --nodes=10
#SBATCH --exclusive          #Number of compute nodes to allocate
###SBATCH --account=m.2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=2017-12-49     #Hans

#load the gromacs module
export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 320 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2016.5-plumed_2.3.5

###REST2###
#First run
aprun $APRUN_OPTIONS gmx_mpi mdrun -plumed plumed.dat -s topol -deffnm md_r -multi 10 -replex 1000 -hrex -maxh 23.5

##Continue Run
# aprun $APRUN_OPTIONS gmx_mpi mdrun -plumed plumed.dat -s topol -cpi md_r -deffnm md_r -multi 10 -replex 1000 -hrex -maxh 23.5
