import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

with open('results.dat', 'r') as f:
    n,x,y = np.loadtxt(f, delimiter='\t', comments='#', usecols=(0,1,2), unpack=True, \
          dtype={'names':('Name', 'DDG', 'Resid_time'), \
                 'formats':('S8', np.float, np.float)})

fig, ax = plt.subplots()

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
line = slope * x + intercept

ax.plot(x,y,'bo')
ax.plot(x,line,'r-',linewidth=2)

ax.set(xlabel=r'$\Delta \Delta G\ (kcal/mol)$', ylabel='Residence time (ns)')

for i,txt in enumerate(n):
    ax.annotate(txt,(x[i],y[i]))

fig.savefig('ASEM_Analogues_Results.png')
