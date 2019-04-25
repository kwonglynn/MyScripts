# USAGE
# python compare.py

# import the necessary packages
from skimage.measure import compare_ssim as ssim
import matplotlib.pyplot as plt
import numpy as np
import cv2
import glob

def calc_ssim(imageA, imageB, color):
    if color == 'rgb':
        sim_image = ssim(imageA, imageB, multichannel = True)
    elif color == 'grey':
        imageA = cv2.cvtColor(imageA, cv2.COLOR_BGR2GRAY)
        imageB = cv2.cvtColor(imageB, cv2.COLOR_BGR2GRAY)
        sim_image = ssim(imageA, imageB, multichannel = False)
    
    return sim_image

## Number of images files:
N = len(glob.glob('fes-*.png')) 

sim_list = []
for i in range(N - 1):
    imageA = 'fes-{}.png'.format(i)
    imageB = 'fes-{}.png'.format(i + 1)
    imageA = cv2.imread(imageA)
    imageB = cv2.imread(imageB)
    sim_image = calc_ssim(imageA, imageB, color = 'rgb')
    sim_list.append(sim_image)

times = []
for i in range(1, N):
    times.append(i)

times = np.array(times) / 10

fo = open('FES_similarity.csv', 'w')
fo.write('Time(ms),Similarity\n')
for i in range(0, N - 1):
    fo.write("{},{.2f}\n".format(times[i],sim_list[i]))

fo.close()
    
fig = plt.figure(figsize = (12,6), dpi = 600)
ax = fig.add_subplot(1, 1, 1)
ax.plot(times, sim_list, 'bo', times, sim_list, 'k')
fig.savefig("FES_similarity_with_time.png")
