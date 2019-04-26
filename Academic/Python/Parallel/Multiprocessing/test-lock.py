# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 10:16:54 2019

@author: Guanglin Kuang
"""

from multiprocessing import Process, Lock
import time

def f1(i, l):
    l.acquire()
    try:
        print('hello world', i)
        time.sleep(1)
    finally:
        l.release()

def f2(i):
    print('hello world', i)
    time.sleep(1)


if __name__ == '__main__':
    lock = Lock()

    for num in range(10):
        Process(target=f1, args=(num,lock)).start()
        #Process(target=f2, args=(num,)).start()