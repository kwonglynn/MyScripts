# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:05:01 2019

@author: Guanglin Kuang
"""

import multiprocessing as mp
import glob
import os

fes_list = glob.glob("fes-*.dat")
def calc_FES(fes):
    print("Processing" + fes)
    N = fes.split('.')[0].split('-')[1]
    # Note: Add the path of python to FES_META.py. And make the script excutable: chmod +x.
    os.system("./FES_METAD.py fes-{}.dat -o fes-{}.png".format(N, N))
    
if __name__ == '__main__':
    pool = mp.Pool(32)
    pool.map(calc_FES, fes_list)
    pool.close()
    print("All fes files processed.")