# -*- coding: utf-8 -*-
from multiprocessing import Pool

def square(x):  
    # calculate the square of the value of x
    return x*x

if __name__ == '__main__':

    # Define the dataset
    dataset = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    # Output the dataset
    print ('Dataset: ' + str(dataset))

    # Run this with a pool of 5 agents having a chunksize of 3 until finished
    agents = 5
    chunksize = 1
    with Pool(processes=agents) as pool:
        result = pool.map(square, dataset, chunksize)

    # Output the result
    print ('Result:  ' + str(result))












