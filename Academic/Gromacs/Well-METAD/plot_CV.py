import matplotlib
import matplotlib.pyplot as plt
import numpy as np

with open('HILLS', 'r') as f:
    t,d,h = np.loadtxt(f, delimiter='\t', comments='#', usecols=(0,1,5), unpack=True, \
          dtype={'names':('time', 'distance', 'height'), \
                 'formats':(np.float, np.float, np.float)})

fig1, ax1 = plt.subplots()

ax1.plot(t,d,'bo')
ax1.plot(t,d,'r-',linewidth=1)
ax1.set(xlabel='Time (ns)', ylabel='distance (ns)')

fig1.savefig('Metad_distance.png')

fig2, ax2 = plt.subplots()

ax2.plot(t,h,'bo')
ax2.plot(t,h,'r-',linewidth=1)
ax2.set(xlabel='Time (ns)', ylabel='height (kJ/mol)')

fig2.savefig('Metad_HILLS.png')
