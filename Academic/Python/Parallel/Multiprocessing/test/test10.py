# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:57:58 2019

@author: Guanglin Kuang
"""

import time
import multiprocessing 

def basic_func(x):
    if x == 0:
        return 'zero'
    elif x%2 == 0:
        return 'even'
    else:
        return 'odd'

def multiprocessing_func(x):
    y = x*x
    print('{} squared results in a/an {} number'.format(x, basic_func(y)))
    time.sleep(2)
    
if __name__ == '__main__':
    
    starttime = time.time()
    pool = multiprocessing.Pool(8)
    pool.map(multiprocessing_func, range(0,10))
    pool.close()
    print('That took {} seconds'.format(time.time() - starttime))