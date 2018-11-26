#!/bin/bash

## Name of the script
#SBATCH -J AMBER

## Names the error and output files according to the jobid
#SBATCH --error=AMBER_%J.err
#SBATCH --output=AMBER_%J.out
#SBATCH --mail-type=ALL

#SBATCH --time=3-00:00:00                #Requested walltime. 10h in this case.
#SBATCH --nodes=2 --exclusive          #Number of compute nodes to allocate
##SBATCH --account=snic2013-1-236       
##SBATCH --account=snic2017-12-49

module load Amber/16-nsc1-intel-2018a-eb 

echo Started at `date`

# suit your needs. This is a generic example:

#(1) Minimization with strong restraints.
pmemd.MPI.run -O -i 01_Min.in -o 01_Min.out -p complex_wi.prmtop  -c complex_wi.inpcrd -ref complex_wi.inpcrd -r 01_Min.rst

#(2) Minimization with no restraints.
pmemd.MPI.run -O -i 02_Min.in -o 02_Min.out -p complex_wi.prmtop  -c 01_Min.rst -r 02_Min.rst

#(3) Heating with weak restraints.
pmemd.MPI.run -O -i 03_Heat.in -o 03_Heat.out -p complex_wi.prmtop -c 02_Min.rst -ref 02_Min.rst -r 03_Heat.rst -x 03_Heat.nc

#(4) Density Equilibration with weak restraints.
pmemd.MPI.run -O -i 04_Density.in -o 04_Density.out -p complex_wi.prmtop -c 03_Heat.rst -ref 03_Heat.rst -r 04_Density.rst -x 04_Density.nc

#(5) Production with no restraints.
pmemd.MPI.run -O -i 05_Prod.in -o 05_Prod.out -p complex_wi.prmtop -c 04_Density.rst -ref 04_Density.rst -r 05_Prod.rst -x 05_Prod.nc

#(6) MMGBSA Calculation.
mpprun MMPBSA.py.MPI -O -i 06_MMPBSA.in -o MMPBSA.dat -sp complex_wi.prmtop -cp complex.prmtop -rp protein.prmtop -lp LIG.prmtop -y 05_Prod.nc >& 06_MMPBSA.log



echo Ended at: `date`
exit 0
#END OF SCRIPT
