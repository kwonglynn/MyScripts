# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 16:28:40 2019

@author: Guanglin Kuang
"""

import multiprocessing as mp
import numpy as np
import time

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

def howmany_within_range_rowonly(row, minimum=4, maximum=8):
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

if __name__ == '__main__':
    # Generate the dataset
    np.random.RandomState(100)
    arr = np.random.randint(0, 10, size=(20000000,5))
    data = arr.tolist()
    
    '''
    # Run without parallelization
    results = []
    start_time = time.time()
    
    for row in data:
        results.append(howmany_within_range(row, minimum=4, maximum=8))
        
    end_time = time.time()
    print("Serial time: {:.2f} seconds".format(end_time - start_time))
    '''
    
    '''    
    # Parallel with pool.apply
    start_time = time.time()
    
    pool = mp.Pool(mp.cpu_count())
    results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]
    pool.close()
    
    end_time = time.time()
    
    print("Parallel time: {:.2f} seconds".format(end_time - start_time))
    '''
    
    '''
    # Parallel with pool.map
    pool = mp.Pool(mp.cpu_count())
    start_time = time.time()    
    results = pool.map(howmany_within_range_rowonly, [row for row in data])
    pool.close()    
    end_time = time.time()    
    print("Parallel time: {:.2f} seconds".format(end_time - start_time))
    '''
    
    '''
    # Parallel with pool.map, in the with mode
    agents = 8
    chunksize = 1
    start_time = time.time()
    with mp.Pool(processes=agents) as pool:
        result = pool.map(howmany_within_range_rowonly, [row for row in data], chunksize)

    end_time = time.time()
    
    print("Parallel time: {:.2f} seconds".format(end_time - start_time))
    '''
    
    #'''
    # Parallel with pool.map_async
    pool = mp.Pool(8)
    start_time = time.time() 
    pool.map(howmany_within_range_rowonly, data)
    pool.close()
    pool.join()
    end_time = time.time()
    print("Parallel time: {:.2f} seconds".format(end_time - start_time))
    #'''
    
    
    