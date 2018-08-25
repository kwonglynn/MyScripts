import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('HILLS', \
                   sep='\s*', engine='python', \
                   comment='#', header=None, \
                   names=['time', 'dist', 'dihed', 'sigma_dist', 'sigma_dihed', 'height', 'biasf'])

time = np.array(data['time']) / 1000
dist = np.array(data['dist']) * 10
height = np.array(data['height']) / 4.184

fig1, ax1 = plt.subplots()

ax1.plot(time, dist, 'bo')
ax1.set(xlabel='Time (ns)', ylabel=r'$distance\ (\AA)$')

fig1.savefig('Metad_distance.png')

fig2, ax2 = plt.subplots()
ax2.plot(time, height, 'bo')
ax2.set(xlabel='Time (ns)', ylabel='height (kcal/mol)')

fig2.savefig('Metad_HILLS.png')
