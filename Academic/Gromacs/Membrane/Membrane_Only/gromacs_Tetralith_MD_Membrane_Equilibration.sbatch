#!/bin/bash -l

## Name of the script
#SBATCH -J MD

## Names the error and output files according to the jobid
#SBATCH --error=MD_%J.err
#SBATCH --output=MD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=2-00:00:00                #Requested walltime 
#SBATCH --nodes=10
#SBATCH --exclusive          #Number of compute nodes to allocate
####SBATCH --account=snic2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=snic2017-12-49     #Hans

#load the gromacs module
module load GROMACS/2018.1-PLUMED-nsc2-gcc-2018a-eb
run=`which gmx_mpi`

###EM###
mkdir EM
$run grompp -f em_membrane.mdp -c POPS_wi.gro -r POPS_wi.gro -p topol.top -o EM/em.tpr
cd EM
mpprun $run mdrun -deffnm em
cd ..

###NVT###
mkdir NVT
$run grompp -f nvt_membrane.mdp -c EM/em.gro -r EM/em.gro -p topol.top -n index.ndx -o NVT/nvt.tpr
cd NVT
mpprun $run mdrun -deffnm nvt
cd ..

###NPT1###
mkdir NPT1
$run grompp -f npt1_membrane.mdp -c NVT/nvt.gro -r NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT1/npt1.tpr
cd NPT1
mpprun $run mdrun -deffnm npt1
cd ..

###NPT2###
mkdir NPT2
$run grompp -f npt2_membrane.mdp -c NPT1/npt1.gro -r NPT1/npt1.gro -t NPT1/npt1.cpt -p topol.top -n index.ndx -o NPT2/npt2.tpr
cd NPT2
mpprun $run mdrun -deffnm npt2
cd ..

###MD###
mkdir MD
$run grompp -f md_membrane.mdp -c NPT2/npt2.gro -t NPT2/npt2.cpt -p topol.top -n index.ndx -o MD/md1.tpr
cd MD
mpprun $run mdrun -deffnm md1
cd ..

