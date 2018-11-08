#!/bin/bash -l

## Name of the script
#SBATCH -J PATH

## Names the error and output files according to the jobid
#SBATCH --error=METAD_%J.err
#SBATCH --output=METAD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=23:59:50                #Requested walltime 
#SBATCH --ntasks 256
#SBATCH --nodes=8
#SBATCH --exclusive          #Number of compute nodes to allocate
###SBATCH --account=m.2017-1-425       #Tu, Account string for the project that you wish to account the job to
#SBATCH --account=2017-12-49     #Hans

#load the gromacs module
export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 256 -d 1 -cc none"
module load cdt/17.10
module swap PrgEnv-cray PrgEnv-gnu
module load plumed/2.5b
module load gromacs/2018.3-plumed_2.5b

###METAD###
##First Run
# gmx_seq grompp -f metad_protein.mdp -c md1.gro -t md1.cpt -p topol.top -n index.ndx -o METAD/metad.tpr
# aprun $APRUN_OPTIONS gmx_mpi mdrun -deffnm metad -plumed plumed_3CV_PATH_FAMetad0.dat

##Continue Run
aprun $APRUN_OPTIONS gmx_mpi mdrun -s metad.tpr -cpi metad.cpt -deffnm metad -plumed plumed_3CV_PATH_FAMetadR.dat
