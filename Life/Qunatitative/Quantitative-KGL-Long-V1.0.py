# -*- coding: utf-8 -*-
# Notice: Change time before usage!
# Note with Tstart, cannot get the latest data!

"""
Created on Sun Aug 27 00:39:55 2017
@author: Guanglin Kuang

#Reference1: http://tushare.org/classifying.html
#Reference2: http://www.waditu.cn/trading.html#id3
"""

import tushare as ts
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']

##### First Part:
##### Definition of the functions.

#################Definition for simple moving average (SMA)#################################
def CalSMA(Close,k):
    SMA = pd.Series(0.0, index = Close.index)
    for i in range(k-1,len(Close)):
        SMA[i] = sum(Close[(i-k+1):(i+1)]/k)
    return SMA

def EWMACal(Close, period, expo):
    EWMA = pd.Series(0.0, index = Close.index)
    EWMA[period-1] = np.mean(Close[:period])
    for i in range(period, len(Close)):
        EWMA[i] = expo*Close[i]+(1-expo)*EWMA[i-1]
    return EWMA   

def GetStock(code, Tperiod):
    while True:
        if TEST == 'True':
            stock = ts.get_k_data(code, ktype = Tperiod, start=Tstart, end=Tend)
        else:
            stock = ts.get_k_data(code, ktype= Tperiod, start=Tstart)
        
        try:
            if len(stock) > 0:
                return stock
        except NameError:
            continue  

def ProperStock(stock_info, code, Tperiod):
    # Pre-filter with basic data.
    if all([float(stock_info.loc[stock_info.index == code].pe) < 200, \
        float(stock_info.loc[stock_info.index == code].pb) < 20, \
        float(stock_info.loc[stock_info.index == code].perundp) > 0]):
#        float(stock_info.loc[stock_info.index==code].rev) < -20, \
#        float(stock_info.loc[stock_info.index==code].profit) < -20]):
        # float(stock_info.loc[stock_info.index==code].gpr) < 15, \
        # float(stock_info.loc[stock_info.index==code].npr) < 0, ]  
        
        # The stock should have enough data for the monthly analysis, > 50 months. 
        stock = GetStock(code, Tperiod)
        if len(stock) > 1050:
            return True
        else:
            return None
    # Exit if the basic contiditons are not met.
    else:
        return None

#MACD is most important!
def MACD(code, Tperiod): 
    # Get the date depending the time period, namely weekly or monthly.
    stock = GetStock(code, Tperiod)

    stock.index = stock.iloc[:,0]
    stock.index = pd.to_datetime(stock.index, format = '%Y-%m-%d')
    stock = stock.iloc[:,1:-1]
    Close = stock.close
            
    DIFF = EWMACal(Close,12,2.0/(1+12)) - EWMACal(Close,26,2.0/(1+26))
    DEA = EWMACal(DIFF,9,2.0/(1+9))
    
    return DIFF, DEA

def KDJ(code, Tperiod):
    #Weekly KDJ.
    stock = GetStock(code, Tperiod)

    stock.index = stock.iloc[:,0]
    stock.index = pd.to_datetime(stock.index, format = '%Y-%m-%d')
    stock = stock.iloc[:,1:-1]
    Close = stock.close
    High = stock.high
    Low = stock.low
    
    dates = Close.index.to_series()    
    periodHigh = pd.Series(np.zeros(len(Close)), index = Close.index)
    periodLow = pd.Series(np.zeros(len(Close)), index = Close.index)
    RSV = pd.Series(np.zeros(len(Close)), index = Close.index)

    for i in range(len(Close)-50, len(Close)):
        period = dates[i-8:i+1]
        periodHigh[i] = High[period].max()
        periodLow[i] = Low[period].min()
        RSV[i] = 100*(Close[i]-periodLow[i])/(periodHigh[i]-periodLow[i])
    
    KValue = pd.Series(50.0,index = Close.index)
    for i in range(len(Close)-50,len(Close)):
        KValue[i] = 2.0/3*KValue[i-1]+1.0/3*RSV[i]
        
    DValue = pd.Series(50.0, index = RSV.index)
    for i in range(len(Close)-50, len(RSV)):
        DValue[i] = 2.0/3*DValue[i-1] + 1.0/3*KValue[i]

    return KValue, DValue
        
##### Second Part:
##### Process all the stocks one by one
if __name__ == '__main__':
    today = datetime.date.today().strftime('%Y-%m-%d')
    fo = open('Stocks-%s.txt' % today,'w')
    #The start time has to be specified explicitly, otherwise the default data is not enough for monthly analysis.
    Tstart = '2012-01-01'

    ######################TEST#########################
    TEST = 'False'       
    if TEST == 'True':
        code_test = '300059'
    #    Tstart=''
    #    Tend=''
        Tstart = '2012-01-01'
        Tend = '2018-05-21'
    ####################TEST##########################
    
    #Get a list of all the stock codes.
    while True:
        stock_info = ts.get_stock_basics()   
        try:
            if len(stock_info) > 0:
                break
        except NameError:
            continue
    
    codes = sorted(list(stock_info.index))
 
    ######################TEST#########################
    if TEST == 'True':
        codes = []
        codes.append(code_test)
     ######################TEST#########################
   
    M = 0
    codes.append('300059')
    print (codes[-10:])      
    for code in codes[-10:]:   
        M += 1
        print(M)
        
        #Pre-Filter
        if not ProperStock(stock_info, code, 'D'):
            continue
        print ('Here')
        # Monthly KDJ and MACD
        Tperiod = 'M'
        KValue, DValue = KDJ(code, Tperiod)
        if not all([KValue[-1] >= DValue[-1], DValue[-1] >= DValue[-2]]):
            continue
        
        DIFF, DEA = MACD(code, Tperiod)
        if not all([DIFF[-1] >= DEA[-1], DEA[-1] >= DEA[-2]]):
            continue
        
        # Weekly KDJ and MACD:
        Tperiod = 'W'
        KValue, DValue = KDJ(code, Tperiod)
        if not all([KValue[-1] >= DValue[-1], DValue[-1] >= DValue[-2]]):
            continue
        
        DIFF, DEA = MACD(code, Tperiod)
        if not all([DIFF[-1] >= DEA[-1], DEA[-1] >= DEA[-2]]):
            continue    
        
        print("%s, %s\n" %(code, today))
        fo.write("%s, %s\n" %(code, today))
       
    fo.close()         
