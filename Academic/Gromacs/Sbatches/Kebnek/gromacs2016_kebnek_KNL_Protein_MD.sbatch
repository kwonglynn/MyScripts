#!/bin/bash 

## Name of the script
#SBATCH -J Keb_KNL 
## Names the error and output files according to the jobid
#SBATCH --error=EMNVTNPTMD_%J.err
#SBATCH --output=EMNVTNPTMD_%J.out
#SBATCH --time=2-00:00:00                #Requested walltime 
###SBATCH --account=SNIC2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=SNIC2017-12-49       #Hans
# Ask for  tasks for use by MPI
#SBATCH -n 280
# Request 4 hardware threads per core
#SBATCH --threads-per-core=4
# Ask for core and memory layout
#SBATCH --constraint=hemi,cache
# Ask for knl queue
#SBATCH -p knl

# It is always best to do a ml purge before loading modules in a submit file
ml purge > /dev/null 2>&1

# Load the module for GROMACS and its prerequisites.
#module load GCC/5.4.0-2.26 CUDA/8.0.61_375.26 impi/2017.3.196
#module load GROMACS/2016.3-PLUMED

ml GCC/5.4.0-2.26 CUDA/8.0.61_375.26 impi/2017.3.196
ml GROMACS/2016.3-PLUMED

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

# Automatic selection of ntomp argument based on "-c" argument to sbatch.
# GROMACS on our KNLs usually runs fastest using only 2 hardware threads per core
# The recommendation is therefore to NOT specify any -c argument to sbatch
# But just in case it is done this code handles that too
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    ntomp="-ntomp $SLURM_CPUS_PER_TASK"
    OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    # Allocate two hardware threads per core and pin them.
    ntomp="-pin on -pinoffset 0 -pinstride 2 -ntomp 2"
    OMP_NUM_THREADS=2
fi
export OMP_NUM_THREADS

###EM###
mkdir EM
$GMX grompp -f em_protein.mdp -c 5KXI_2FA_wi.gro -p topol.top -o EM/em.tpr
cd EM
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm em
cd ..

###NVT###
mkdir NVT
$GMX grompp -f nvt_protein.mdp -c EM/em.gro -p topol.top -n index.ndx -o NVT/nvt.tpr
cd NVT
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm nvt
cd ..

###NPT###
mkdir NPT
$GMX grompp -f npt_protein.mdp -c NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT/npt.tpr
cd NPT
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm npt
cd ..

###MD###
mkdir MD
$GMX grompp -f md_protein.mdp -c NPT/npt.gro -t NPT/npt.cpt -p topol.top -n index.ndx -o MD/md1.tpr
cd MD
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm md1
cd ..

