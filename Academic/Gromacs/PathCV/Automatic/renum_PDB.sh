for PDB in `cat chainA-loopC.txt`; do
    name=$(basename $PDB .pdb)
    gmx_mpi editconf -f $PDB -o ${name}_renum.pdb -resnr 185
done
