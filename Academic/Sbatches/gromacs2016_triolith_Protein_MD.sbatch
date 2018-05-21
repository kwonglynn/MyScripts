#!/bin/bash -l

## Name of the script
#SBATCH -J 5KXI_2FA

## Names the error and output files according to the jobid
#SBATCH --error=EMNVTNPTMD_%J.err
#SBATCH --output=EMNVTNPTMD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=2-00:00:00                #Requested walltime 
#SBATCH --nodes=8 
#SBATCH --exclusive          #Number of compute nodes to allocate
####SBATCH --account=snic2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=snic2017-12-49     #Hans

#load the gromacs module
module add gromacs/5.1.4-plumed
run=`which gmx_mpi`

###EM###
mkdir EM
$run grompp -f em_protein.mdp -c 5KXI_2FA_wi.gro -p topol.top -o EM/em.tpr
cd EM 
mpprun $run mdrun -deffnm em
cd ..

###NVT###
mkdir NVT
$run grompp -f nvt_protein.mdp -c EM/em.gro -p topol.top -n index.ndx -o NVT/nvt.tpr 
cd NVT
mpprun $run mdrun -deffnm nvt
cd ..

###NPT###
mkdir NPT
$run grompp -f npt_protein.mdp -c NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT/npt.tpr 
cd NPT
mpprun $run mdrun -deffnm npt
cd ..

###MD###
mkdir MD
$run grompp -f md_protein.mdp -c NPT/npt.gro -t NPT/npt.cpt -p topol.top -n index.ndx -o MD/md1.tpr
cd MD
mpprun $run mdrun -deffnm md1
cd ..
