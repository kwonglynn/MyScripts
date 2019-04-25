# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 17:20:22 2019

@author: Guanglin Kuang
"""

import os
import time
from multiprocessing import Process
 
def doubler(number):
    """
    A doubling function that can be used by a process
    """
    result = number * 2
    proc = os.getpid()
    print('{0} doubled to {1} by process id: {2}'.format(number, result, proc))
    time.sleep(5)
 
if __name__ == '__main__':
    numbers = [5, 10, 15, 20, 25]
    procs = []
 
    for index, number in enumerate(numbers):
        proc = Process(target=doubler, args=(number,))
        procs.append(proc)
        proc.start()
 
    for proc in procs:
        proc.join()