'''
Reference:
1, https://www.shanelynn.ie/python-pandas-read_csv-load-data-from-csv-files/
2, https://chrisalbon.com/python/data_wrangling/pandas_dataframe_importing_csv/  
3, https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html
4, https://matplotlib.org/tutorials/introductory/pyplot.html#sphx-glr-tutorials-introductory-pyplot-py
5, https://stepik.org/lesson/44680/step/1?unit=22723
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'

data = pd.read_csv('HILLS', \
                   sep='\s*', engine='python', \
                   comment='#', header=None, \
                   names=['time', 'dist', 'dihed', 'sigma_dist', 'sigma_dihed', 'height', 'biasf'])

time = np.array(data['time']) / 1000
dist = np.array(data['dist']) * 10
height = np.array(data['height']) / 4.184

fig1, ax1 = plt.subplots()
ax1.plot(time, dist, 'k-', linewidth=1)
plt.xlim(0, 5000)
plt.ylim(0, 70)
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel(r'$distance\ (\AA)$', fontsize=12)
plt.tick_params(labelsize=10)
fig1.savefig('Metad_distance.png', dpi=600)

fig2, ax2 = plt.subplots()
ax2.plot(time, height, 'k-', linewidth=1)
plt.xlim(0, 5000)
plt.ylim(0.00, 0.14)
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel('height (kcal/mol)', fontsize=12)
plt.tick_params(labelsize=10)
fig2.savefig('Metad_HILLS.png', dpi=600)