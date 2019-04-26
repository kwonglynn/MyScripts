# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:28:37 2019

@author: Guanglin Kuang
"""

from multiprocessing import Pool
import time 

def doubler(number):
    return number * 2
 
if __name__ == '__main__':
    numbers = range(50000000)
    start_time = time.time()
    
    pool = Pool(processes=8)
    result = pool.map(doubler, numbers)
    
    end_time = time.time()
    
    print("Parallel time: {:.2f} seconds".format(end_time - start_time))