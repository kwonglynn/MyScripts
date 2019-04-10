# Load the refere structure, which should be loaded first.
mol new MIT1_NMR.pdb
# use this structure for the reference
set reference [atomselect top backbone]

# Load the trajectory
mol new npt_remd0.gro
mol addfile md_remd0.xtc waitfor all

set file [open "rmsd.dat" w]

set compare [atomselect 1 backbone]
set all [atomselect 1 all]
set num_frames [molinfo 1 get numframes]
for {set frame 0} {$frame < $num_frames} {incr frame} {
   # get the correct frame
   $compare frame $frame
   # compute the transformation
   set trans_mat [measure fit $compare $reference]
   # do the alignment
   $all frame $frame
   $all move $trans_mat
   # compute the RMSD
   set rmsd [measure rmsd $compare $reference]
   # print the RMSD
   set t [ expr {$frame/100.0}]
   puts -nonewline $file "[format %-6.3f $t]\t"
   puts $file "[format %6.2f $rmsd]"
}  
flush $file
close $file
