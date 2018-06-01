# !/usr/bin/python
# -*- coding: utf-8 -*-
from urllib.request import urlopen
import tushare as ts
import numpy as np
from scipy import stats
import datetime

##################################Main#######################################
today = datetime.date.today().strftime('%Y-%m-%d')
fo = open('Fund_Mining-%s.txt' % today, 'w')
fo.write("Code\tYear1\tYear2\tYear3\tSlope\tSharp\tSTD\n")

StartDate = '2014-05-31'
EndDate = '2018-05-31'

##Extract the equity funds, dump monetary and bond funds.
ResponseAllFunds = urlopen('http://fund.eastmoney.com/js/fundcode_search.js')
AllFundsTxt = str(ResponseAllFunds.read())

#处理数据 将其转化为list
AllFundsTxt = AllFundsTxt[AllFundsTxt.find('=')+2:AllFundsTxt.rfind(';')]
AllFundsList = eval(AllFundsTxt)

#AllFundsList = [['110022']]       # For TEST only:

# 循环处理每个基金, 只保留需要的基金。
i = 0
for fund in AllFundsList:
    i += 1
#    if i > 5: break            # For TEST only!
    code = str(fund[0])
    print (code)
    try:
        fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
        prices = np.array(fund.total.dropna())
        returns = np.array(fund.change.dropna())
        # 1. New funds. 2. Monetary funds. 3. Bond funds. 4. Bad funds. 5. Not good funds.
        if any([len(prices) < 751, \
               fund.value[0] > fund.total[0], \
               fund.total[0] < 1.1, \
               prices[0] < prices[-1], \
               (prices[0]-prices[250])/prices[250] < 0.2]):
            continue
    except:
        continue
    print ('Found:', code)
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


    
    
    

