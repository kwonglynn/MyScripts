# !/usr/bin/python
# -*- coding: utf-8 -*-
import tushare as ts
import numpy as np
from scipy import stats
import datetime

##################################Main########################################
today = datetime.date.today().strftime('%Y-%m-%d')
fo = open('Fund_Mining-%s.txt' % today, 'w')
fo.write("Code\tYear1\tYear2\tYear3\tSlope\tSharp\tSTD\n")

StartDate = '2014-05-31'
EndDate = '2018-05-31'
 
TEST = 'False'
if TEST == 'True': 
    codes = ['110022','540006','590008','519606','090013','020026']

if TEST != 'True':
    with open('EquityFunds-%s.dat' % today,'r') as fi:
        codes=fi.readlines()
    
for code in codes:
    code = code.strip()
    print (code)
    try:
        fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
        prices = np.array(fund.total.dropna())
        returns = np.array(fund.change.dropna())
    
        if len(prices) < 751:
            continue
        elif prices[0] < prices[-1] or (prices[0]-prices[250])/prices[250] < 0.2:
            continue
    except:
        continue

    #Return in one and two years
    year1 = (prices[0]-prices[250])/prices[250]*100
    year2 = (prices[250]-prices[500])/prices[500]*100
    year3 = (prices[500]-prices[750])/prices[750]*100
    
    x = np.arange(len(prices))
    slope = -1 * stats.linregress(x, prices)[0] * 100
    
    #SHARP
    ER = np.mean(returns)
    Rf = 5.0 / 250
    std = np.std(returns)
    sharp = (ER-Rf)/std * 100

    fo.write("'%s'\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n" %  (code, year1, year2, year3, slope, sharp, std))
    
fo.close()
        

    
    
    

