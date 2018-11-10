import numpy as np

MSDs = []
with open('rmsd.dat', 'r') as fi:
    for line in fi:
        MSD = np.square(float(line.strip()))
        MSDs.append(MSD)

mean_MSD = np.mean(MSDs)
LAMBDA = 2.3/mean_MSD

print (LAMBDA)
