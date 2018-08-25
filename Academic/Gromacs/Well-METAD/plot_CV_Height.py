import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('HILLS', \
                   sep='\s*', engine='python', \
                   comment='#', header=None, \
                   names=['time', 'dist', 'dihed', 'sigma_dist', 'sigma_dihed', 'height', 'biasf'])

time = np.array(data['time'])
dist = np.array(data['dist'])
height = np.array(data['height'])

fig1, ax1 = plt.subplots()

ax1.plot(time, dist, 'bo')
ax1.set(xlabel='Time (ns)', ylabel='distance (ns)')

fig1.savefig('Metad_distance.png')

fig2, ax2 = plt.subplots()
ax2.plot(time, height, 'bo')
ax2.set(xlabel='Time (ns)', ylabel='height (kJ/mol)')

fig2.savefig('Metad_HILLS.png')
