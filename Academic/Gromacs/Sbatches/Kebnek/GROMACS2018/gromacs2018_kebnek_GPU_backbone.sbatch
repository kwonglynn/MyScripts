#!/bin/bash
# Change to your actual SNIC project number
#SBATCH -A SNIC2018-2-38
# Name of the job
#SBATCH -J MD-GPU
# Asking for walltime
#SBATCH -t 2-00:00:00
# Use for 4 tasks
#SBATCH -n 4
# Usefor 7 threads per task
#SBATCH -c 7
# Remember the total number of cores = 28 (kebnekaise) so that
# n x c = 28 (running on a single node, or multiples of 28 for multi node runs)
# Asking for 2 K80 GPU cards per node
#SBATCH --gres=gpu:k80:2

# It is always best to do a ml purge before loading modules in a submit file
ml purge > /dev/null 2>&1

# Load the module for GROMACS and its prerequisites.
# Use spider to check if the GOMACS version is GPU-enabled.
module load GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
module load GROMACS/2018.3-PLUMED

# Automatic selection of single or multi node based GROMACS
GMX="gmx_mpi"
MPIRUN="mpirun"
ntmpi=""

# Automatic selection of ntomp argument based on "-c" argument to sbatch
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    ntomp="$SLURM_CPUS_PER_TASK"
else
    ntomp="1"
fi
# Make sure to set OMP_NUM_THREADS equal to the value used for ntomp
# to avoid complaints from GROMACS
export OMP_NUM_THREADS=$ntomp

###EM###
mkdir EM
gmx grompp -f em_protein.mdp -c complex_wi.gro -r complex_wi.gro -p complex.top -o EM/em.tpr
cd EM
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm em
cd ..

###NVT###
mkdir NVT
gmx grompp -f nvt_protein.mdp -c EM/em.gro -r EM/em.gro -p complex.top -n index.ndx -o NVT/nvt.tpr
cd NVT
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm nvt
cd ..

###NPT###
mkdir NPT
gmx grompp -f npt_protein.mdp -c NVT/nvt.gro -r NVT/nvt.gro -t NVT/nvt.cpt -p complex.top -n index.ndx -o NPT/npt.tpr
cd NPT
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm npt
cd ..

###MD###
mkdir MD
cp plumed-MD.dat MD
gmx grompp -f md_backbone.mdp -c NPT/npt.gro -r NPT/npt.gro -t NPT/npt.cpt -p complex.top -n index.ndx -o MD/md.tpr
cd MD
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm md -plumed plumed-MD.dat
cd ..

