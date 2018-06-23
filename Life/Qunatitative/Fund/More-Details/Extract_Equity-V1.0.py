# !/usr/bin/python
# -*- coding: utf-8 -*-
#import datetime
import urllib2
import json
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
#import re
import tushare as ts

def ExtractFund(i, codes):
    ResponseAllFunds = urllib2.urlopen('http://fund.eastmoney.com/js/fundcode_search.js')
    AllFundsTxt = ResponseAllFunds.read()

    #处理数据 将其转化为list
    AllFundsTxt = AllFundsTxt[AllFundsTxt.find('=')+2:AllFundsTxt.rfind(';')]
    AllFundsList = json.loads(AllFundsTxt.decode('utf-8'))

    # 循环处理每个基金, 只保留需要的基金。
    with open('AllFunds.dat','w') as fo:
        for fund in AllFundsList:
            i += 1
            print 'i: ', i
#            if i > 10: break
            code = str(fund[0])
            fo.write(code+'\n')
            codes.append(code)

def FilterFund(j, StartDate, EndDate, codes):
    with open('EquityFunds.dat','w') as fo:
        for code in codes:
            print code
            try:
                fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
                if float(fund.value) <= float(fund.total):
                    j += 1
                    print 'j: ', j                
                    fo.write(code+'\n')
            except AttributeError:
                continue
                
###############################User Input############################################
StartDate = '2018-01-19'
EndDate = '2018-01-19'
i=0
codes = []
ExtractFund(i, codes)
FilterFund(i, StartDate, EndDate, codes)

