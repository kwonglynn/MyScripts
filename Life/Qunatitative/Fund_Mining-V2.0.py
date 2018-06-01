# !/usr/bin/python
# -*- coding: utf-8 -*-
import tushare as ts
import numpy as np
from scipy import stats

##################################Main########################################
StartDate = '2013-05-31'
EndDate = '2018-05-31'
 
TEST = 'True'
if TEST == 'True': 
    codes = ['110022','540006','590008','519606','090013','020026']

fo = open('Fund_Mining2.dat','w')
fo.write("Code\tSlope\tSTD\n")

if TEST != 'True':
    with open('EquityFunds-old.dat','r') as fi:
        codes=fi.readlines()
    
for code in codes:
    code = code.strip()
    fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
    prices = np.array(fund.total.dropna())
    returns = np.array(fund.change.dropna())
    
    if len(prices) < 500:
        continue
    
    print (code)
    x = np.arange(len(prices))
    slope = stats.linregress(x, prices)[0] * 100
    std = np.std(returns)

    fo.write("%s\t%4.2f\t%4.2f\n" %  (code, slope, std))
    
fo.close()
        

    
    
    

