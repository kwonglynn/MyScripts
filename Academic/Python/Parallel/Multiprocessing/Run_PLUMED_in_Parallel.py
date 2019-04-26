# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:05:01 2019

@author: Guanglin Kuang
"""

import multiprocessing as mp
import glob
import os

def run_plumed(hills):
    print("Processing: " + hills)
    N = hills.split('-')[1]
    os.system("plumed sum_hills --hills HILLS-{} --bin 98,99 --outfile fes-{}.dat".format(N, N))
    
if __name__ == '__main__':
    hills_list = glob.glob("HILLS-*")
    pool = mp.Pool(32)
    pool.map(run_plumed, hills_list)
    pool.close()
    print("All HILLS files processed.")