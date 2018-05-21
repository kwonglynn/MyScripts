# -*- coding: utf-8 -*-
# Reference:
# https://pythonspot.com/matplotlib-histogram/
# https://matplotlib.org/examples/statistics/histogram_demo_features.html

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab

fi1 = open('EPP_r0.xvg','r')
fi2 = open('EPP_r9.xvg','r')

all = []
for i in range(10):
    pots = []
    fi = open("EPP_r%d.xvg" % i)
    for line in fi:
        if not line.startswith('#') and not line.startswith('@'):
            pot = float(line.strip().split()[1])/4.184
            pots.append(pot)
            all.append(pot)
  
    Apots = np.array(pots)
      
    mean = np.mean(Apots)
    variance = np.var(Apots)
    sigma = np.sqrt(variance)
    x = np.linspace(min(Apots), max(Apots), 200)
#    n, bins, patches = plt.hist(Apots, 100, density=True, histtype='step', alpha=0.75)
#    n, bins, patches = plt.hist(Apots, 100,  alpha=0.75, normed=True, cumulative=True)

    n, bins, patches = plt.hist(Apots, 200, normed=True, histtype='step', alpha=0.75)
    plt.plot(x, mlab.normpdf(x, mean, sigma), label="%s" % i)
    fi.close()
  
plt.xlim((min(all), max(all)))
plt.xlabel('Potential Energy (kcal/mol)')
#plt.xlabel(r'$E_{PP}\ (kcal/mol)$')
plt.ylabel('Probability')
plt.legend()
plt.show()