from molecule import *
from AtomSel import AtomSel

load('complex_box.gro')
CA =  AtomSel('name CA')

print CA

quit()
