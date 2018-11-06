#!/bin/bash -l

## Name of the script
#SBATCH -J MD

## Names the error and output files according to the jobid
#SBATCH --error=MD_%J.err
#SBATCH --output=MD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=7-00:00:00                #Requested walltime 
#SBATCH --nodes=10
#SBATCH --exclusive          #Number of compute nodes to allocate
####SBATCH --account=snic2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=snic2017-12-49     #Hans

#load the gromacs module
module load GROMACS/2018.1-PLUMED-nsc2-gcc-2018a-eb
run=`which gmx_mpi`

###EM###
mkdir EM
$run grompp -f em_membrane.mdp -c complex_wi.gro -r complex_wi.gro -p topol.top -o EM/em.tpr
cd EM
mpprun $run mdrun -deffnm em
cd ..

###NVT###
mkdir NVT
$run grompp -f nvt_membrane.mdp -c EM/em.gro -r EM/em.gro -p topol.top -n index.ndx -o NVT/nvt.tpr
cd NVT
mpprun $run mdrun -deffnm nvt
cd ..

###NPT###
mkdir NPT
$run grompp -f npt_membrane.mdp -c NVT/nvt.gro -r NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT/npt.tpr
cd NPT
mpprun $run mdrun -deffnm npt
cd ..

###MD###
mkdir MD
$run grompp -f md_membrane.mdp -c NPT/npt.gro -t NPT/npt.cpt -p topol.top -n index.ndx -o MD/md1.tpr
cd MD
mpprun $run mdrun -deffnm md1
cd ..

