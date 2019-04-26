# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:13:36 2019

@author: Guanglin Kuang
"""

from multiprocessing import Process, Lock
import time
import os
 
 
def printer(item, lock):
    """
    Prints out the item that was passed in
    """
    lock.acquire()
    proc = os.getpid()

    try:
        print(item)
        print('Process id: {}'.format( proc))
        time.sleep(5)
    finally:
        lock.release()
    

 
if __name__ == '__main__':
    lock = Lock()
    items = range(11)
    for item in items:
        p = Process(target=printer, args=(item, lock))
        p.start()