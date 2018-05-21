antechamber -fi pdb -i LIG.pdb -fo mol2 -o LIG_bcc_gaff.mol2 -c bcc -at gaff -rn LIG -nc 1
parmchk -f mol2 -i LIG_bcc_gaff.mol2 -o LIG.frcmod
# Process with tleap
tleap -s -f leap_gaff.in
