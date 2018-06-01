# !/usr/bin/python
# -*- coding: utf-8 -*-
#import datetime
from urllib.request import urlopen
import tushare as ts
import datetime
#import json
#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')
#import re
         
###############################User Input############################################
StartDate = '2018-05-31'
EndDate = '2018-05-31'
today = datetime.date.today().strftime('%Y-%m-%d')

ResponseAllFunds = urlopen('http://fund.eastmoney.com/js/fundcode_search.js')
AllFundsTxt = str(ResponseAllFunds.read())

#处理数据 将其转化为list
AllFundsTxt = AllFundsTxt[AllFundsTxt.find('=')+2:AllFundsTxt.rfind(';')]
AllFundsList = eval(AllFundsTxt)

# 循环处理每个基金, 只保留需要的基金。
i = 0
codes = []
with open('AllFunds-%s.dat' % today,'w') as fo:
    for fund in AllFundsList:
        i += 1
        print ('i: ', i)
#            if i > 10: break
        code = str(fund[0])
        fo.write(code + '\n')
        codes.append(code)

j = 0
with open('EquityFunds-%s.dat' % today,'w') as fo:
    for code in codes:
        print (code)
        try:
            fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
            if float(fund.value) <= float(fund.total) or (int(fund.total)*10) > 11:
                j += 1
                print ('j: ', j)
                fo.write(code + '\n')
        except:
            continue

