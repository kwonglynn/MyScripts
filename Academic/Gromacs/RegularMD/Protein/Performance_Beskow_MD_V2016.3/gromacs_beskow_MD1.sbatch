#!/bin/bash -l

## Name of the script
#SBATCH -J MD

## Names the error and output files according to the jobid
#SBATCH --error=EMNVTNPTMD_%J.err
#SBATCH --output=EMNVTNPTMD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=00:10:00                #Requested walltime 
#SBATCH --ntasks YYY
#SBATCH --nodes=XXX 
#SBATCH --exclusive          #Number of compute nodes to allocate
#SBATCH --account=m.2017-1-425       #Tu, Account string for the project that you wish to account the job to
###SBATCH --account=2017-12-49     #Hans

#load the gromacs module
export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n YYY -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2016.3-plumed_2.3

###MD Continue###
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm md1
