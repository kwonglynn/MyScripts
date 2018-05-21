#!/bin/bash -l

# Batch script to run a GPU job on G02 with SGE queuing system.

# 1. Set the name of the job (Only change here for amber GPU jobs).
#$ -N CD73_MD

# 2. Force bash as the executing shell.
#$ -S /bin/bash

# 3. Request a number of GPU cards:
#$ -l gpus=1

# 4. Request the wallclock time (format hours:minutes:seconds).
#$ -l h_rt=100:00:0

# 5. Merge the standard error stream into the standard output stream.
#$ -j y

# 6. Request a queue for the job.
#$ -q gpu.q

# 7. Make sure the current environment variables are used on the SGE jobs.
#$ -V

# 8. Use the current working directory for input and output.
#$ -cwd 

echo Started at: `date`

export CUDA_VISIBLE_DEVICES=1,2

# 9. Run the Amber Jobs:
#################################################################Job Section############################################################
#(1) Minimization with strong restraints.
pmemd.cuda -O -i 01_Min.in -o 01_Min.out -p complex_wi.prmtop  -c complex_wi.inpcrd -r 01_Min.rst -ref complex_wi.inpcrd

#(2) Minimization with no restraints.
pmemd.cuda -O -i 02_Min.in -o 02_Min.out -p complex_wi.prmtop  -c 01_Min.rst -r 02_Min.rst

#(3) Heating with weak restraints.
mpirun -np 2 pmemd.cuda.MPI -O -i 03_Heat.in -o 03_Heat.out -p complex_wi.prmtop -c 02_Min.rst -r 03_Heat.rst -ref 02_Min.rst -x 03_Heat.nc

#(4) Density Equilibration with weak restraints.
mpirun -np 2 pmemd.cuda.MPI -O -i 04_Density.in -o 04_Density.out -p complex_wi.prmtop -c 03_Heat.rst -r 04_Density.rst -ref 03_Heat.rst -x 04_Density.nc

#(5) Production with no restraints.
mpirun -np 2 pmemd.cuda.MPI -O -i 05_Prod.in -o 05_Prod.out -p complex_wi.prmtop -c 04_Density.rst -r 05_Prod.rst -x 05_Prod.nc
########################################################################################################################################

# 10. Submit the job with this command:
# qsub amber_SGE_G02.sh

# 11. Check the status of the job with qstat or qmon (GUI).

echo Ended at: `date`
exit 0
#END OF SCRIPT
