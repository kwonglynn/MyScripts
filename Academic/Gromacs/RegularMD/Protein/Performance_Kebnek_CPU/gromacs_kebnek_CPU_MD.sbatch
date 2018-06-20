#!/bin/bash 

## Name of the script
#SBATCH -J Keb_CPU 
## Names the error and output files according to the jobid
#SBATCH --error=EMNVTNPTMD_%J.err
#SBATCH --output=EMNVTNPTMD_%J.out
#SBATCH --time=00:10:00                #Requested walltime 
###SBATCH --account=SNIC2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=SNIC2017-12-49       #Hans
# Ask for  tasks for use by MPI
#SBATCH -n YYY
# Ask for 2 cores per MPI-task
#SBATCH -c 2

# It is always best to do a ml purge before loading modules in a submit file
ml purge > /dev/null 2>&1

# Load the module for GROMACS and its prerequisites.
module load GCC/5.4.0-2.26 CUDA/8.0.61_375.26 impi/2017.3.196
module load GROMACS/2016.3-PLUMED

# Automatic selection of single or multi node based GROMACS
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then
    GMX="gmx_mpi"
    MPIRUN="mpirun"
    ntmpi=""
else
    GMX="gmx"
    MPIRUN=""
    ntmpi="-ntmpi $SLURM_NTASKS"
fi

# Automatic selection of ntomp argument based on "-c" argument to sbatch
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    ntomp="$SLURM_CPUS_PER_TASK"
else
    ntomp="1"
fi
# Make sure to set OMP_NUM_THREADS equal to the value used for ntomp
# to avoid complaints from GROMACS
export OMP_NUM_THREADS=$ntomp

###MD###
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm md1

