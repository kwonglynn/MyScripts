#!/bin/bash -l

## Name of the script
#SBATCH -J METAD

## Names the error and output files according to the jobid
#SBATCH --error=METAD_%J.err
#SBATCH --output=METAD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=7-00:00:00                #Requested walltime 
#SBATCH --nodes=10
#SBATCH --exclusive          #Number of compute nodes to allocate
####SBATCH --account=snic2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=snic2017-12-49     #Hans

#load the gromacs module
module load GROMACS/2018.1-PLUMED-nsc2-gcc-2018a-eb
run=`which gmx_mpi`

##METAD Run
# $run grompp -f metad_backbone.mdp -c NPT/npt.gro -t NPT/npt.cpt -p complex.top -n index.ndx -o METAD/metad1.tpr
#mpprun $run mdrun -deffnm metad1 -plumed plumed.dat -maxh 167

##Continue Run
mpprun $run mdrun -s metad1.tpr -cpi metad1.cpt -deffnm metad1 -plumed plumed.dat -maxh 167


