#!/bin/bash -l

## Name of the script
#SBATCH -J Trio_MD

## Names the error and output files according to the jobid
#SBATCH --error=REST2_%J.err
#SBATCH --output=REST2_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=7-00:00:00                #Requested walltime 
#SBATCH --nodes=10 
#SBATCH --exclusive          #Number of compute nodes to allocate
#SBATCH --account=snic2017-1-425       #Tu, Account string for the project that you wish to account the job to
###SBATCH --account=snic2017-12-49     #Hans

#load the gromacs module
module load gromacs/2016.5-plumed2.3.5
run=`which gmx_mpi`

###REST2###
# $run grompp -f rest2_protein.mdp -c NPT/npt.gro -t NPT/npt.cpt -p topol.top -n index.ndx -o topol.tpr
mpprun $run mdrun -plumed plumed.dat -s topol -deffnm md_r -multi 10 -replex 1000 -hrex -maxh 167
