# Load the refere structure, which should be loaded first.
mol new MIT1_NMR.pdb
# use this structure for the reference
set reference [atomselect top backbone]

for {set n 0} {$n < 31} {incr n} {
   #Preprocess 
   echo 1 | trjconv_mpi -s md_remd${n}.tpr -f md_remd${n}.xtc -pbc mol -ur compact -o md_remd${n}.noPBC.noWater.xtc

   # Load the trajectory
   mol new MIT1_processed.gro
   mol addfile md_remd${n}.noPBC.noWater.xtc waitfor all

   set file [open "REMD${n}_rmsd.dat" w]

   set compare [atomselect top backbone]
   set all [atomselect top all]

   set num_frames [molinfo top get numframes]
   for {set frame 0} {$frame < $num_frames} {incr frame} {
      # get the correct frame
      $compare frame $frame

      # compute the transformation
      set trans_mat [measure fit $compare $reference]
      # do the alignment
      set all [atomselect top all]
      $all frame $frame
      $all move $trans_mat

      # compute the RMSD
      set rmsd [measure rmsd $compare $reference]
      # print the RMSD
      set t [ expr {$frame/1000.0}]
      puts -nonewline $file "[format %-6.2f $t]\t"
      puts $file "[format %6.2f $rmsd]"
    }  
    flush $file
    close $file
}
quit
