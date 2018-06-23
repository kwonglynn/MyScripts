# -*- coding: utf-8 -*-
# Notice: Change time before usage!
# Note with Tstart, cannot get the latest data!
TEST='False'
    
if TEST=='True':
    code_test='601318'
#    Tstart=''
#    Tend=''
    Tstart='2015-11-16'
    Tend='2017-04-28'

"""
Created on Sun Aug 27 00:39:55 2017
@author: Guanglin Kuang
"""
#Reference1: http://tushare.org/classifying.html
#Reference2: http://www.waditu.cn/trading.html#id3

import tushare as ts
import pandas as pd
import numpy as np
import datetime
import time
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']
ts.set_token('e39db8588179dd3d9744dcc6925e83c978068f6617cb19e0b52abc349afb4a7c')

##### First Part:
##### Definition of the functions.

#################Definition for simple moving average (SMA)#################################
def CalSMA(Close,k):
    SMA=pd.Series(0.0, index=Close.index)
    for i in range(k-1,len(Close)):
        SMA[i]=sum(Close[(i-k+1):(i+1)]/k)
    return SMA

def EWMACal(Close, period, expo):
    EWMA=pd.Series(0.0, index=Close.index)
    EWMA[period-1]=np.mean(Close[:period])
    for i in range(period, len(Close)):
        EWMA[i]=expo*Close[i]+(1-expo)*EWMA[i-1]
    return EWMA

###############################Definition for Wu Lian Yang##################################
def WuLianYang(Open, Close, High, Low, Signal):
    for i in range(len(Close)-10, len(Close)-5):
        Yang=[]
        incr=[]
        Kincrease1=(Close[i+1]-Close[i])/Close[i]
        Kvibrate1=(High[i+1]-Low[i+1])/Low[i+1]
        Kextend1=(Close[i+1]-Open[i+1])/Open[i+1]
        if all([Kincrease1 >= 0, Kincrease1 < 0.03, Kvibrate1 < 0.05, Kextend1 > -0.005]):
            Yang.append(1)
            
        Kincrease2=(Close[i+2]-Close[i+1])/Close[i+1]
        Kvibrate2=(High[i+2]-Low[i+2])/Low[i+2]
        Kextend2=(Close[i+2]-Open[i+2])/Open[i+2]
        if all([Kincrease2 >= 0, Kincrease2 < 0.03, Kvibrate2 < 0.05, Kextend2 > -0.005]):
            Yang.append(1)
        if Low[i+2] >= Low[i+1]:
            incr.append(1)
            
        Kincrease3=(Close[i+3]-Close[i+2])/Close[i+2]
        Kvibrate3=(High[i+3]-Low[i+3])/Low[i+3]
        Kextend3=(Close[i+3]-Open[i+3])/Open[i+3]
        if all([Kincrease3 >= 0, Kincrease3 < 0.03, Kvibrate3 < 0.05, Kextend3 > -0.005]):
            Yang.append(1)
        if Low[i+3] >= Low[i+2]:
            incr.append(1)
            
        Kincrease4=(Close[i+4]-Close[i+3])/Close[i+3]
        Kvibrate4=(High[i+4]-Low[i+4])/Low[i+4]
        Kextend4=(Close[i+4]-Open[i+4])/Open[i+4]
        if all([Kincrease4 >= 0, Kincrease4 < 0.03, Kvibrate4 < 0.05, Kextend4 > -0.005]):
            Yang.append(1)
        if Low[i+4] >= Low[i+3]:
            incr.append(1)
        
        Kincrease5=(Close[i+5]-Close[i+4])/Close[i+4]
        Kvibrate5=(High[i+5]-Low[i+5])/Low[i+5]
        Kextend5=(Close[i+5]-Open[i+5])/Open[i+5]
        if all([Kincrease5 >= 0, Kincrease5 < 0.03, Kvibrate5 < 0.05, Kextend5 > -0.005]):
            Yang.append(1)
        if Low[i+5] >= Low[i+4]:
            incr.append(1)
            
        if all([len(Yang) >= 3, len(incr) >=2, (Close[i+5]-Close[i+1])/Close[i+1] > 0.01, \
                Kincrease1 > -0.02, Kincrease1 < 0.03, Kvibrate1 < 0.06, Kextend1 > -0.02, Kextend1 < 0.03, \
                Kincrease2 > -0.02, Kincrease2 < 0.03, Kvibrate2 < 0.06, Kextend2 > -0.02, Kextend2 < 0.03, \
                Kincrease3 > -0.02, Kincrease3 < 0.03, Kvibrate3 < 0.06, Kextend3 > -0.02, Kextend3 < 0.03, \
                Kincrease4 > -0.02, Kincrease4 < 0.03, Kvibrate4 < 0.06, Kextend4 > -0.02, Kextend4 < 0.03, \
                Kincrease5 > -0.02, Kincrease5 < 0.03, Kvibrate5 < 0.06, Kextend5 > -0.02, Kextend5 < 0.03]):
            Signal[i+5]=1            

#####################################Definition for Jin Cha###########################################
def JinCha(Close, SMA20, SMA120, Signal):
    for i in range(len(Close)-10,len(Close)):
        if all([abs(SMA20[i]-SMA120[i])/SMA120[i] < 0.005, \
                SMA20[i] >= SMA20[i-1], 
                SMA120[i] >= SMA120[i-1], \
                Close[i] >= SMA20[i], Close[i] >= SMA120[i]]):
            Signal[i]=1

#######################################Definition for MA Duo Tou#####################################
def DuoTou(Close, SMA5, SMA10, SMA20,SMA60, SMA120, SMA250, Signal):
    for i in range(len(Close)-10,len(Close)-3):
        if all([Close[i] >= SMA5[i], (Close[i]-SMA5[i])/SMA5[i] < 0.1, \
                Close[i] >= SMA10[i], (Close[i]-SMA10[i])/SMA10[i] < 0.12, \
                Close[i] >= SMA20[i], (Close[i]-SMA20[i])/SMA20[i] < 0.15, \
                Close[i] >= SMA60[i], (Close[i]-SMA60[i])/SMA60[i] < 0.2, \
                Close[i] >= SMA120[i], (Close[i]-SMA120[i])/SMA120[i] < 0.25, \
                Close[i] >= SMA250[i], (Close[i]-SMA250[i])/SMA250[i] < 0.3, \
                (Close[i+1]-Close[i])/Close[i] > -0.05, \
                (Close[i+2]-Close[i])/Close[i] > -0.05, \
                (Close[i+3]-Close[i])/Close[i] > -0.05]):
            Signal[i+3]=1

#Important!
#####################################Definition for K Tu Pou#########################################
def TuPou(Close, Volume, Signal):
    for i in range(len(Close)-10,len(Close)-2):
        Kchange1=(Close[i+1]-Close[i])/Close[i]
        Vchange1=(Volume[i+1]-Volume[i])/Volume[i]
        Kchange2=(Close[i+2]-Close[i+1])/Close[i+1]
          
        if all([Kchange1 >= 0.035, Vchange1 >= 0.2 , Kchange2 >= -0.02]):
            Signal[i+2] = 1
                   
####################################Definition for Zhang Ting#######################################
def ZhangTing(Close, Signal):
    for i in range(len(Close)-30,len(Close)-1):     
        if all([(Close[i]-Close[i-1])/Close[i-1] > 0.095, \
                (Close[i+1]-Close[i])/Close[i]  < 0.095]):
            Signal[i+1] = 1

#######################################Duo Fang Xin Hao#####################################
def DuoFangXinHao(Open, Close, High, Low, Signal):
    for i in range(len(Close)-10, len(Close)-3):
        Kchange1=(Close[i+1]-Close[i])/Close[i]
        Kextend1=(Close[i+1]-Open[i+1])/Open[i+1]
        Kvibrate1=(High[i+1]-Low[i+1])/Low[i+1]
        
        Kchange2=(Close[i+2]-Close[i+1])/Close[i+1]
        Kextend2=(Close[i+2]-Open[i+2])/Open[i+2]

        Kchange3=(Close[i+3]-Close[i+2])/Close[i+2]
        Kextend3=(Close[i+3]-Open[i+3])/Open[i+3]
        
        # Duo Fang Pao  
        if all([Kchange1 < -0.03, Kextend1 < -0.03, Kextend2 > 0.03]):
            Signal[i+2] = 1
        
        # Shang Zhang Bao Xian
        if all([Kchange1 < -0.02, Kextend1 > -0.03, Kextend1 < 0.01, Kvibrate1 < 0.05, \
                (Open[i+2]-Low[i+1])/Low[i+1] < -0.01, (Close[i+2]-High[i+1])/High[i+1] > 0.01, \
                Kchange2 > 0.01, Kextend2 > 0.03]):
            Signal[i+2] = 1
        
        # Chui Zi Xian
        com=min(Open[i+2], Close[i+2])
        if all([(Low[i+2]-com)/com < -0.03, \
                Low[i+2] < Low[i+1], \
                Low[i+3] > Low[i+2]]):
            Signal[i+3] = 1

        # Zao Chen Zhi Xin
        if all([Kchange1 < -0.03, Kextend1 < -0.03, \
                Kchange2 > -0.02, Kchange2 < 0.01, abs(Kextend2) < 0.02,\
                Kchange3 > 0.03, Kextend3 > 0.03]):
            Signal[i+3]=1

        # TiaoKong
        if all([(Low[i+3]-High[i+2])/High[i+2] > 0.01, \
                (Low[i+3]-High[i+2])/High[i+2] < 0.99]):
            Signal[i+3] = 1

####################################Definition for Wen Jian#######################################
def WenJian(Close, SMA5, Signal):
    for i in range(len(Close)-10, len(Close)-4):
        Yang=[]
        incr=[]

        if Close[i+1] >= SMA5[i+1]:
            Yang.append(1)
        if Low[i+1] >= SMA5[i+1]:
            incr.append(1)

        if Close[i+2] >= SMA5[i+2]:
            Yang.append(1)
        if Low[i+2] >= SMA5[i+2]:
            incr.append(1)
            
        if Close[i+3] >= SMA5[i+3]:
            Yang.append(1)
        if Low[i+3] >= SMA5[i+3]:
            incr.append(1)

        if Close[i+4] >= SMA5[i+4]:
            Yang.append(1)
        if Low[i+4] >= SMA5[i+4]:
            incr.append(1)            

        if all([len(Yang) >= 3, len(incr) >=2]):
            Signal[i+4]=1        

#MACD is most important!
def MACD(Close, DIFF, DEA, Signal):    
    #Di Wei Jin Cha
    for i in range(len(Close)-10, len(Close)-1):        
        if all([DIFF[i]<= DEA[i], DIFF[i]<1.0, DIFF[i+1]>=DEA[i+1]]):
            Signal[i+1]=1
    
    #Chuan Guo 0 Zhou
    for i in range(len(Close)-10, len(Close)-2):
        if all([DEA[i+1]<=0, DEA[i+2]>=0]):
            Signal[i+2]=1

    #Jiang Cha Wei Cha
    for i in range(len(Close)-10, len(Close)-2):
        if all([DIFF[i]>=DEA[i], DIFF[i+1]<=DIFF[i], DIFF[i+1]<=DIFF[i+2], \
                abs(DIFF[i+1]-DEA[i+1])<=0.1, DIFF[i+2] >= DEA[i+2], \
                DEA[i]<=DEA[i+1], DEA[i+1]<=DEA[i+2]]):
            Signal[i+2]=1       

    # Nianhe
    for i in range(len(Close)-10, len(Close)-5):
        Yang=[]
        incr=[]

        if abs(DIFF[i+1]-DEA[i+1])<=0.1:
            Yang.append(1)
        if (DIFF[i+1]-DEA[i+1])>=0 and (DIFF[i+1]-DEA[i+1])<=0.1:
            incr.append(1)        

        if abs(DIFF[i+2]-DEA[i+2])<=0.1:
            Yang.append(1)
        if (DIFF[i+2]-DEA[i+2])>=0 and (DIFF[i+2]-DEA[i+2])<=0.1:
            incr.append(1)            

        if abs(DIFF[i+3]-DEA[i+3])<=0.1:
            Yang.append(1)            
        if (DIFF[i+3]-DEA[i+3])>=0 and (DIFF[i+3]-DEA[i+3])<=0.1:
            incr.append(1)            

        if abs(DIFF[i+4]-DEA[i+4])<=0.1:
            Yang.append(1)
        if (DIFF[i+4]-DEA[i+4])>=0 and (DIFF[i+4]-DEA[i+4])<=0.1:
            incr.append(1)           

        if abs(DIFF[i+5]-DEA[i+5])<=0.1:
            Yang.append(1) 
        if (DIFF[i+5]-DEA[i+5])>=0 and (DIFF[i+5]-DEA[i+5])<=0.1:
            incr.append(1)           
           
        if all([len(Yang) >=4, len(incr) >=3]):
            Signal[i+5]=1     

def DKDJ(Close):
    #Dayly KDJ.
    dates=Close.index.to_series()
    
    periodHigh=pd.Series(np.zeros(len(Close)), index=Close.index)
    periodLow=pd.Series(np.zeros(len(Close)), index=Close.index)
    RSV=pd.Series(np.zeros(len(Close)), index=Close.index)

    for i in range(len(Close)-50, len(Close)):
        period=dates[i-8:i+1]
        periodHigh[i]=High[period].max()
        periodLow[i]=Low[period].min()
        RSV[i]=100*(Close[i]-periodLow[i])/(periodHigh[i]-periodLow[i])
    
    KValue=pd.Series(50.0,index=Close.index)
    for i in range(len(Close)-50,len(Close)):
        KValue[i]=2.0/3*KValue[i-1]+1.0/3*RSV[i]
        
    DValue=pd.Series(50.0,index=RSV.index)
    for i in range(len(Close)-50,len(RSV)):
        DValue[i]=2.0/3*DValue[i-1]+1.0/3*KValue[i]
    
    if all([DValue[-2] <= DValue[-1], KValue[-1] >= DValue[-1]]):
        return True
    else:
        return False

def WKDJ(code):
    #Weekly KDJ.
    try:
        if TEST=='True':
            stock2=ts.get_k_data(code, ktype = 'W', start=Tstart, end=Tend)
        else:
            stock2=ts.get_k_data(code, ktype = 'W')
    except:
        time.sleep(1)
        if TEST=='True':
            stock2=ts.get_k_data(code, ktype = 'W', start=Tstart, end=Tend)
        else:
            stock2=ts.get_k_data(code, ktype = 'W')

    stock2.index=stock2.iloc[:,0]
    stock2.index=pd.to_datetime(stock2.index,format='%Y-%m-%d')
    stock2=stock2.iloc[:,1:-1]
    Close2=stock2.close
    
    dates2=Close2.index.to_series()    
    periodHigh=pd.Series(np.zeros(len(Close2)), index=Close2.index)
    periodLow=pd.Series(np.zeros(len(Close2)), index=Close2.index)
    RSV=pd.Series(np.zeros(len(Close2)), index=Close2.index)

    for i in range(len(Close2)-50, len(Close2)):
        period=dates2[i-8:i+1]
        periodHigh[i]=High[period].max()
        periodLow[i]=Low[period].min()
        RSV[i]=100*(Close2[i]-periodLow[i])/(periodHigh[i]-periodLow[i])
    
    KValue=pd.Series(50.0,index=Close2.index)
    for i in range(len(Close2)-50,len(Close2)):
        KValue[i]=2.0/3*KValue[i-1]+1.0/3*RSV[i]
        
    DValue=pd.Series(50.0,index=RSV.index)
    for i in range(len(Close2)-50,len(RSV)):
        DValue[i]=2.0/3*DValue[i-1]+1.0/3*KValue[i]

    if all([KValue[-2] <= KValue[-1], DValue[-2] <= DValue[-1], \
            KValue[-1] >= DValue[-1]]):
        return True
    else:
        return False
        
##### Second Part:
##### Process all the stocks one by one
stock_info=ts.get_stock_basics()
codes=sorted(list(stock_info.index))

today=datetime.date.today().strftime('%Y-%m-%d')
fo=open('Stocks-%s.txt' % today,'w')
    
M=0

if TEST=='True':
    codes=[]
    codes.append(code_test)
        
for code in codes:   
    M+=1
    print M
    
    if float(stock_info.loc[stock_info.index==code].pe) > 200:
        continue
    elif float(stock_info.loc[stock_info.index==code].pb) > 15:
        continue
    elif float(stock_info.loc[stock_info.index==code].esp) < 0:
        continue
    elif float(stock_info.loc[stock_info.index==code].perundp) < 0:
        continue
    elif float(stock_info.loc[stock_info.index==code].rev) < -20:
        continue
    elif float(stock_info.loc[stock_info.index==code].profit) < -20:
        continue

#    elif float(stock_info.loc[stock_info.index==code].gpr) < 15:
#        continue
#    elif float(stock_info.loc[stock_info.index==code].npr) < 0:
#        continue        

    try:
        if TEST=='True':
            stock=ts.get_k_data(code, start=Tstart, end=Tend)
        else:
            stock=ts.get_k_data(code)
    except:
        time.sleep(1)
        if TEST=='True':
            stock=ts.get_k_data(code, start=Tstart, end=Tend)
        else:
            stock=ts.get_k_data(code)

    if len(stock) < 250:
        continue

    stock.index=stock.iloc[:,0]
    stock.index=pd.to_datetime(stock.index,format='%Y-%m-%d')
    stock=stock.iloc[:,1:-1]
    Open=stock.open
    Close=stock.close
    High=stock.high
    Low=stock.low
    Volume=stock.volume    

###Main filter 1:
    DIFF=EWMACal(Close,12,2.0/(1+12)) - EWMACal(Close,26,2.0/(1+26))
    DEA=EWMACal(DIFF,9,2.0/(1+9))

    if not all([(round(DIFF[-2],2) <= round(DIFF[-1],2) \
                 or round(DEA[-2],2) <= round(DEA[-1],2)), \
                DEA[-1]>=-0.5, DEA[-1]<=5.0]):
        continue

###Main filter 2:
    SMA20=CalSMA(Close,20)
    SMA250=CalSMA(Close,250)
    if not all([round(Close[-1],2) >= round(SMA20[-1],2), round(SMA20[-2],2) <= round(SMA20[-1],2), \
                round(Close[-1],2) >= round(SMA250[-1],2), round(SMA250[-4],2) <= round(SMA250[-1],2)]):
        continue

###Delicate filters:
    Max=max(Close)
    Min=min(Close)  
    
    SMA5=CalSMA(Close,5)
    SMA10=CalSMA(Close,10)
    SMA60=CalSMA(Close,60)
    SMA120=CalSMA(Close,120)
   
    if not all([(Close[-1]-SMA10[-1])/SMA10[-1] >= -0.05, (Close[-1]-SMA10[-1])/SMA10[-1] <= 0.1, \
            round(SMA10[-2],2) <= round(SMA10[-1],2), \
            round(SMA60[-2],2) <= round(SMA60[-1],2), \
            round(SMA120[-3],2) <= round(SMA120[-1],2), \
            (Close[-1]-SMA20[-1])/SMA20[-1] >= 0, (Close[-1]-SMA20[-1])/SMA20[-1] <= 0.12, \
            (Close[-1]-SMA60[-1])/SMA60[-1] >= 0, (Close[-1]-SMA60[-1])/SMA60[-1] <= 0.22, \
            (max(Close[-10:])-min(Close[-10:]))/min(Close[-60:]) <= 0.2, \
            (max(Close[-60:])-min(Close[-60:]))/min(Close[-60:]) <= 0.6]):
        continue

    Signal=pd.Series(0,index=Close.index)
    WuLianYang(Open, Close, High, Low, Signal)
    JinCha(Close, SMA20, SMA120, Signal)
    DuoTou(Close, SMA5, SMA10, SMA20,SMA60, SMA120, SMA250, Signal)
    TuPou(Close, Volume, Signal)
    ZhangTing(Close, Signal)
    DuoFangXinHao(Open, Close, High, Low, Signal)
    WenJian(Close, SMA5, Signal)
    MACD(Close, DIFF, DEA, Signal)    

#####Proces and write out the results.
    for i in range(len(Close)-10, len(Close)):
        if Signal[i]==1 and DKDJ(Close) and WKDJ(code):
            fo.write("%s\n" %(code))
            break
 
# Write out the signal date, for testing.
    if TEST=='True':
        print code        
        for i in range(len(Close)-10, len(Close)):
            if Signal[i]==1 and DKDJ(Close) and WKDJ(code):
                BuyIn=Signal[i:i+1]
                BuyIn=str(BuyIn).split('\n')[1:-1]
                
                for line in BuyIn:
                    print line         
   
fo.close()         
