# Load the refere structure, which should be loaded first.
mol new MIT1_NMR_DS.pdb

# Select the backbone atoms of the reference structure. The backbone atoms will be used for RMSD calculation in this case.
set sel_ref [atomselect top backbone]

# Read the file listing the target structures.
set file_list [open models.list]
set content [read $file_list]
close $file_list

# Create a new file to write the calculated results.
set file_rmsd [open models.rmsd w]

# Split into records on newlines. Split will return a list in TCL type.
set targets [split $content "\n"]

# Iterative over the records, the target structures should be processed one by one.
foreach target $targets {
	# This is to deal with the last blank list member.
	if { $target == ""} {break}
	# Load a structure to be compared.
	mol new $target

	# Select the backbone atoms of the target structure.
	set sel_fit [atomselect top backbone]

	# Select all the atoms of the target structure. All the atoms should be moved in order to align, not just the backbone atoms.
	set sel_all [atomselect top all]

	# Compute the transformation matrix between the target structure and reference structure.
	set M [measure fit $sel_fit $sel_ref]
	$sel_all move $M

	# Compute the RMSD.
	set rmsd [measure rmsd $sel_fit $sel_ref]

	# Write the calculated RMSD of each target structure to a result file.
	puts $file_rmsd "$target \t $rmsd"
	
	# Unload the current structure before the next iteration.
	mol delete top
}
close $file_rmsd
