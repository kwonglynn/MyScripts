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

fig1, ax1 = plt.subplots(111)
ax1.plot(time, dist, 'k-', linewidth=1)
plt.xlim(0, 5000)
plt.ylim(0, 70)
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel(r'$distance\ (\AA)$', fontsize=12)
plt.tick_params(labelsize=10)
fig1.savefig('Metad_distance.png', dpi=600)

fig2, ax2 = plt.subplots(111)
ax2.plot(time, height, 'k-', linewidth=1)
plt.xlim(0, 5000)
plt.ylim(0.00, 0.14)
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel('height (kcal/mol)', fontsize=12)
plt.tick_params(labelsize=10)
fig2.savefig('Metad_HILLS.png', dpi=600)