#!/bin/bash
# Change to your actual SNIC project number
#SBATCH -A SNIC2017-12-49
# Asking for walltime
#SBATCH -J MD-CPU
#SBATCH -t 7-00:00:00
# Use 14 tasks
#SBATCH -n 140
# Use 2 threads per task
#SBATCH -c 2

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
gmx grompp -f em_membrane.mdp -c POPS_wi.gro -r POPS_wi.gro -p topol.top -o EM/em.tpr
cd EM
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm em
cd ..

###NVT###
mkdir NVT
gmx grompp -f nvt_membrane.mdp -c EM/em.gro -r EM/em.gro -p topol.top -n index.ndx -o NVT/nvt.tpr
cd NVT
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm nvt
cd ..

###NPT1###
mkdir NPT1
gmx grompp -f npt1_membrane.mdp -c NVT/nvt.gro -r NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT1/npt1.tpr
cd NPT1
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm npt1
cd ..

###NPT2###
mkdir NPT2
gmx grompp -f npt2_membrane.mdp -c NPT1/npt1.gro -r NPT1/npt1.gro -t NPT1/npt1.cpt -p topol.top -n index.ndx -o NPT2/npt2.tpr
cd NPT2
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm npt2
cd ..

###MD###
mkdir MD
gmx grompp -f md_membrane.mdp -c NPT2/npt2.gro -t NPT2/npt2.cpt -p topol.top -n index.ndx -o MD/md1.tpr
cd MD
$MPIRUN $GMX mdrun $ntmpi -ntomp $ntomp -deffnm md1
cd ..

