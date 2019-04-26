# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:06:33 2019

@author: Guanglin Kuang
"""

import os
 
from multiprocessing import Process, current_process
 
 
def doubler(number=2):
    """
    A doubling function that can be used by a process
    """
    result = number * 2
    proc_name = current_process().name
    print('{0} doubled to {1} by: {2}'.format(number, result, proc_name))
 
 
if __name__ == '__main__':
    numbers = [5, 10, 15, 20, 25]
    procs = []
     
    for index, number in enumerate(numbers):
        proc = Process(target=doubler, args=(number,))
        procs.append(proc)
        proc.start()
 
    proc = Process(target=doubler, name='Test', args=(2,))
    proc.start()
    procs.append(proc)
 
    for proc in procs:
        proc.join()