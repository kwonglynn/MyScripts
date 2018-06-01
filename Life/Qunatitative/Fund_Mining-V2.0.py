# !/usr/bin/python
# -*- coding: utf-8 -*-
import tushare as ts
import numpy as np
from scipy import stats

##################################Main########################################
StartDate = '2015-05-31'
EndDate = '2018-05-31'
 
TEST = 'True'
if TEST == 'True': 
    codes = ['110022','540006','590008','519606','090013','020026']

fo = open('Fund_Mining2.dat','w')
fo.write("Code\tyear1\tyear2\tSlope\tSharp\tSTD\n")

if TEST != 'True':
    with open('EquityFunds-old.dat','r') as fi:
        codes=fi.readlines()
    
for code in codes:
    code = code.strip()
    print (code)
    try:
        fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
    except AttributeError:
        continue
    except ValueError:
        continue
    
    prices = np.array(fund.total.dropna())
    returns = np.array(fund.change.dropna())
    
    if len(prices) < 500:
        continue
    elif prices[0] < prices[-1] or (prices[0]-prices[250])/prices[250] < 0.2:
        continue

    #Return in one and two years
    year1 = (prices[0]-prices[250])/prices[250]*100
    year2 = (prices[0]-prices[500])/prices[500]*100
    
    x = np.arange(len(prices))
    slope = -1 * stats.linregress(x, prices)[0] * 100
    
    #SHARP
    ER = np.mean(returns)
    Rf = 5.0 / 250
    std = np.std(returns)
    sharp = (ER-Rf)/std * 100

    fo.write("%s\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n" %  (code, year1, year2, slope, sharp, std))
    
fo.close()
        

    
    
    

