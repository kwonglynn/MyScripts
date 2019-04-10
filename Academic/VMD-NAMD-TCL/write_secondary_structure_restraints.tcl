# vmd -dispdev text -e write_secondary_structure_restraints.tcl
# Load and align the strucutre

mol new SV2C_4JRA_processed.gro
set file [open "posre_beta_backbone.itp" w]

puts $file "\[ position_restraints \]"
puts $file "; atom  type      fx      fy      fz"

set beta_backbone [atomselect top "sheet and backbone"]

set N [$beta_backbone num]
set serials [$beta_backbone get serial]

for {set i 0} {$i < $N} {incr i} {
    set serial [lindex $serials $i]
    puts $file [format "%d\t1\t1000\t1000\t1000" $serial]
}

flush $file
close $file

quit