#For Amber
module load amber/16-nsc1
antechamber -fi pdb -i LIG.pdb -fo mol2 -o LIG_bcc_gaff.mol2 -c bcc -at gaff -rn LIG 
parmchk -f mol2 -i LIG_bcc_gaff.mol2 -o LIG.frcmod
# Process with tleap
tleap -s -f leap_gaff.in
acpype -p LIG.prmtop -x LIG.inpcrd
sed -e "1,\$s/1  LIG/1LIG  /g" LIG_GMX.gro > LIG_MD.gro
cp LIG_GMX.top LIG_MD.gro ..
