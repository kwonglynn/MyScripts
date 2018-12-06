#!/bin/bash -l

## Name of the script
#SBATCH -J METAD

## Names the error and output files according to the jobid
#SBATCH --error=METAD_%J.err
#SBATCH --output=METAD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=23:59:50                #Requested walltime 
#SBATCH --ntasks=320
#SBATCH --nodes=10
#SBATCH --exclusive          #Number of compute nodes to allocate
#SBATCH --account=2018-3-549
###SBATCH --account=2017-12-49     #Hans

#load the gromacs module
export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 320 -d 1 -cc none"
module load cdt/17.10
module swap PrgEnv-cray PrgEnv-gnu
module load plumed/2.5b
module load gromacs/2018.3-plumed_2.5b

###METAD###
##First Run
# aprun -n 1 gmx_mpi grompp -f metad_backbone.mdp -c MD/md.gro -r MD/md.gro -t MD/md.cpt -p complex.top -n index.ndx -o METAD1/metad1.tpr
aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm metad1 -plumed plumed.dat -maxh 23.5

##Continue Run
# aprun $APRUN_OPTIONS gmx_mpi mdrun -s metad1.tpr -cpi metad1.cpt -deffnm metad1 -plumed plumed.dat -maxh 23.5
