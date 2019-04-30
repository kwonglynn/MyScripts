# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:45:35 2019

@author: Guanglin Kuang
"""

# import the necessary packages
from skimage.measure import compare_ssim as ssim
import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted
import cv2
import glob
import os
import ray

@ray.remote
def calc_ssim(i, imageA, imageB, color = 'rgb'):
    os.chdir('/proj/molmat/users/x_guaku/Tau/Single/MK6240/METAD_Small_VIP/convergence/ray') 
    imageA = cv2.imread(imageA)
    imageB = cv2.imread(imageB)

    if color == 'rgb':
        ssim_image = ssim(imageA, imageB, multichannel = True)
    elif color == 'grey':
        imageA = cv2.cvtColor(imageA, cv2.COLOR_BGR2GRAY)
        imageB = cv2.cvtColor(imageB, cv2.COLOR_BGR2GRAY)
        ssim_image = ssim(imageA, imageB, multichannel = False)

    return (i, ssim_image)

if __name__ == '__main__':
    # Start ray on one node (the driver):
    # ray start --head --redis-port=6379 --num-cpus=16
    # Start ray on other nodes (the workers):
    # ray start --redis-address 10.24.2.48:6379 --num-cpus=16
    # Stop the ray on a node:
    # ray stop    

    ray.init(redis_address="10.24.2.48:6379")
    ## Number of images files:
    fes_pngs = glob.glob('fes-*.png')
    fes_pngs = natsorted(fes_pngs)
    N = len(fes_pngs)
    
    ## Use the last image as reference.
    ref = 'fes-{}.png'.format(N-1)
    results_id = []
    for i in range(len(fes_pngs)):
        png = fes_pngs[i]
        result = calc_ssim.remote(i, png, ref, 'rgb')
        results_id.append(result)

    results = ray.get(results_id)
    print (results)
    results.sort(key=lambda x: x[0])
    ssim_list = [r for i, r in results]
    print (ssim_list)
    print("All FES images files processed.")

    times = []
    for i in range(1, N+1):
        times.append(i)
    
    times = np.array(times) / 10
    
    fo = open('FES_similarity.csv', 'w')
    fo.write("#Time(ms)\tSimilarity\n")
    for i in range(0, N):
        fo.write("{:.1f}\t{:.2f}\n".format(times[i],ssim_list[i]))
    
    fo.close()
    
    fig = plt.figure(figsize = (12,6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(0, 10)
    ax.set_ylim(0.5, 1.0)
    ax.plot(times, ssim_list, 'bo', times, ssim_list, 'k')
    fig.savefig("FES_similarity_with_time.png", dpi=600)
