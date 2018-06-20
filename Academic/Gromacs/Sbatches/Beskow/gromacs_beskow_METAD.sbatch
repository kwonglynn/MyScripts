#!/bin/bash -l

## Name of the script
#SBATCH -J METAD

## Names the error and output files according to the jobid
#SBATCH --error=EMNVTNPTMD_%J.err
#SBATCH --output=EMNVTNPTMD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=23:59:50                #Requested walltime 
#SBATCH --ntasks 128
#SBATCH --nodes=4
#SBATCH --exclusive          #Number of compute nodes to allocate
###SBATCH --account=m.2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=2017-12-49     #Hans

#load the gromacs module
export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 128 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2016.5-plumed_2.3.5 

###METAD###
##First Run
# aprun -n 1 gmx_mpi grompp -f metad_backbone.mdp -c NPT/npt.gro -t NPT/npt.cpt -p complex.top -n index.ndx -o METAD1/metad1.tpr
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm metad1 -plumed plumed.dat -maxh 23.5

##Continue Run
# aprun $APRUN_OPTIONS gmx_mpi mdrun -s metad1.tpr -cpi metad1.cpt -deffnm metad1 -plumed plumed.dat -maxh 23.5
