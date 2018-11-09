#!/bin/bash -l

## Name of the script
#SBATCH -J MD

## Names the error and output files according to the jobid
#SBATCH --error=MD_%J.err
#SBATCH --output=MD_%J.out
####SBATCH --mail-type=ALL

#SBATCH --time=4-00:00:00                #Requested walltime 
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

###NPT1###
mkdir NPT1
$run grompp -f npt1_membrane.mdp -c NVT/nvt.gro -r NVT/nvt.gro -t NVT/nvt.cpt -p topol.top -n index.ndx -o NPT1/npt1.tpr
cp index.ndx protein.pdb plumed_NPT.dat NPT1
cd NPT1
mpprun $run mdrun -deffnm npt1 -plumed plumed_NPT.dat
cd ..

###NPT2###
mkdir NPT2
$run grompp -f npt2_membrane.mdp -c NPT1/npt1.gro -r NPT1/npt1.gro -t NPT1/npt1.cpt -p topol.top -n index.ndx -o NPT2/npt2.tpr
cp index.ndx protein.pdb plumed_NPT.dat NPT2
cd NPT2
mpprun $run mdrun -deffnm npt2 -plumed plumed_NPT.dat
cd ..

###NPT3###
mkdir NPT3
$run grompp -f npt3_membrane.mdp -c NPT2/npt2.gro -r NPT2/npt2.gro -t NPT2/npt2.cpt -p topol.top -n index.ndx -o NPT3/npt3.tpr
cp index.ndx protein.pdb plumed_NPT.dat NPT3
cd NPT3
mpprun $run mdrun -deffnm npt3 -plumed plumed_NPT.dat
cd ..

###MD1###
mkdir MD1
$run grompp -f md1_membrane.mdp -c NPT3/npt3.gro -t NPT3/npt3.cpt -p topol.top -n index.ndx -o MD1/md1.tpr
cp index.ndx protein.pdb plumed_MD1.dat MD1
cd MD1
mpprun $run mdrun -deffnm md1 -plumed plumed_MD1.dat
cd ..

###MD2###
mkdir MD2
$run grompp -f md2_membrane.mdp -c MD1/md1.gro -t MD1/md1.cpt -p topol.top -n index.ndx -o MD2/md2.tpr
cp index.ndx protein.pdb plumed_MD2.dat MD2
cd MD2
mpprun $run mdrun -deffnm md2 -plumed plumed_MD2.dat
cd ..


