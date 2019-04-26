# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:45:35 2019

@author: Guanglin Kuang
"""

# import the necessary packages
from skimage.measure import compare_ssim as ssim
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import cv2
import glob

def calc_ssim(i, imageA, imageB, color = 'rgb'):
    imageA = cv2.imread(imageA)
    imageB = cv2.imread(imageB)

    if color == 'rgb':
        ssim_image = ssim(imageA, imageB, multichannel = True)
    elif color == 'grey':
        imageA = cv2.cvtColor(imageA, cv2.COLOR_BGR2GRAY)
        imageB = cv2.cvtColor(imageB, cv2.COLOR_BGR2GRAY)
        ssim_image = ssim(imageA, imageB, multichannel = False)

    return (i, ssim_image)

def collect_result(result):
    global results
    results.append(result)

if __name__ == '__main__':
    ## Number of images files:
    fes_pngs = glob.glob('fes-*.png')
    N = len(fes_pngs)
    
    ## Use the last image as reference.
    ref = 'fes-{}.png'.format(N-1)
    
    pool = mp.Pool(32)
    results = []
    for i, image in enumerate(fes_pngs):
        pool.apply_async(calc_ssim, args=(image, ref, 'rgb'), callback=collect_result)
    pool.close()
    pool.join()
    
    results.sort(key=lambda x: x[0])
    ssim_list = [r for i, r in results]

    print("All FES images files processed.")

    times = []
    for i in range(1, N+1):
        times.append(i)
    
    times = np.array(times) / 10
    
    fo = open('FES_similarity.csv', 'w')
    fo.write('Time(ms),Similarity\n')
    for i in range(0, N - 1):
        fo.write("{},{}\n".format(times[i],ssim_list[i]))
    
    fo.close()
    
    fig = plt.figure(figsize = (12,6), dpi = 600)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times, ssim_list, 'bo', times, ssim_list, 'k')
    fig.savefig("FES_similarity_with_time.png")
