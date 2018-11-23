echo 0 | gmx trjconv -s md.tpr -f md.gro -o md-noPBC.gro -pbc mol -ur compact -n ../index
echo 24 | gmx trjconv -s md.tpr -f md.xtc -o md-noPBC-noWater.xtc -pbc mol -ur compact -n ../index.ndx
echo -e "5\n2\n" | gmx rms -s md.tpr -f md-noPBC-noWater.xtc -o LIG_rmsd.xvg -tu ns
