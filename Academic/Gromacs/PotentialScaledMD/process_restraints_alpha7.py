num_chain1 = 3298
num_chain2 = 3298

fo_chain1 = open("posre_Protein_scale.itp", 'w')
fo_chain2 = open("posre_Protein2_scale.itp", 'w')

######
import MDAnalysis as mda

u = mda.Universe('complex.gro')

##MDAnalysis requires the topology to have segment IDs
#Need to add segids, otherwise it will be problematic for protein with multiple chains.
chain1 = u.select_atoms("bynum 1:%d or resname LIG" % num_chain1)
chain2 = u.select_atoms("bynum %d:%d or resname LIG" % (num_chain1, num_chain1+num_chain2))

backbone1 = chain1.select_atoms("backbone").ids
free1 = chain1.select_atoms("backbone and (resid 182:192 or (byres (around 6 resname LIG)))").ids
backbone2 = chain2.select_atoms("backbone").ids
free2 = chain2.select_atoms("backbone and (byres (around 6 resname LIG))").ids

print ("Free backbone atoms in the first chain:")
for id in free1: print (*id)
print ("\n")
print ("Free backbone atoms in the second chain:")
for id in free2: print (*id)

title = '''\
[ position_restraints ]
; atom func     g          r            k
'''
fo_chain1.write(title)
fo_chain2.write(title)

func = 2
g = 1
r = 0.1
k = 50
for id in backbone1:
    if id not in free1:
        fo_chain1.write("%d\t%d\t%d\t%f\t%d\n" % (id,func,g,r,k))

for id in backbone2:
    if id not in free2:
        id = id - num_chain1
        fo_chain2.write("%d\t%d\t%d\t%f\t%d\n" % (id,func,g,r,k))

fo_chain1.close()
fo_chain2.close()
