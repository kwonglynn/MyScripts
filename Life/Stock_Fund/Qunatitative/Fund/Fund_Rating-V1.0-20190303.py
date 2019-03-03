# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 19:22:37 2019

@author: Guanglin Kuang
"""

import tushare as ts
import datetime

##################################Main#######################################
today = datetime.date.today().strftime('%Y-%m-%d')
fo = open('Fund_Rating-%s.txt' % today, 'w')
fo.write("Code\tThisYear\tRecent1Y\tRecent2Y\tRecent3Y\tRecent5Y\tY2018\tY2017\tY2016\tY2015\tY2014\n")

# StartDate = '2013-01-01'
# EndDate = '2019-03-03'

AllFundsList = ['519068', '162605', '000628', '000527', '160505', '450002', '020026', '090013', '519069', '519066', \
                '000577', '110011', '180012', '040008', '519732', '000619', '519736', '260108', '540012', '001938', \
                '000742', '000751', '340008', '540006', '070032']

def getValue(date):
    i = 0
    while len(list(prices[date])) == 0:
        i += 1
        date = datetime.datetime.strptime(date, '%Y-%m-%d') - datetime.timedelta(days = i)
        date = date.strftime('%Y-%m-%d')
        
        if datetime.datetime.strptime(date, '%Y-%m-%d') < datetime.datetime.strptime(str(prices.index[-1]).split()[0], '%Y-%m-%d'):
            return prices[-1]
    return prices[date].values[0]

# 循环处理每个基金
for code in AllFundsList:
    print (code)
    fund = ts.fund.nav.get_nav_history(code)
#    fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
    # Use cumulative value
    prices = fund.total.dropna()
    # Use daily value
#    prices = fund.value.dropna()

    #Return in one, two and three years
    ThisYear = (getValue(today) - getValue('2019-01-01')) / getValue('2019-01-01') * 100
    
    ##Yearly return
    OneYearAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365)).strftime('%Y-%m-%d')
    Recent1Y = (getValue(today) - getValue(OneYearAgo)) / getValue(OneYearAgo) * 100

    TwoYearsAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 2)).strftime('%Y-%m-%d')
    Recent2Y = (getValue(today) - getValue(TwoYearsAgo)) / getValue(TwoYearsAgo) * 100

    ThreeYearsAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 3)).strftime('%Y-%m-%d')
    Recent3Y = (getValue(today) - getValue(ThreeYearsAgo)) / getValue(ThreeYearsAgo) * 100

    FiveYearsAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 5)).strftime('%Y-%m-%d')
    Recent5Y = (getValue(today) - getValue(FiveYearsAgo)) / getValue(FiveYearsAgo) * 100

    Y2018 = (getValue('2018-12-31') - getValue('2018-01-01')) / getValue('2018-01-01') * 100
    Y2017 = (getValue('2017-12-31') - getValue('2017-01-01')) / getValue('2017-01-01') * 100
    Y2016 = (getValue('2016-12-31') - getValue('2016-01-01')) / getValue('2016-01-01') * 100
    Y2015 = (getValue('2015-12-31') - getValue('2015-01-01')) / getValue('2015-01-01') * 100
    Y2014 = (getValue('2014-12-31') - getValue('2014-01-01')) / getValue('2014-01-01') * 100


    fo.write("'%s'\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n" % \
             (code,ThisYear,Recent1Y,Recent2Y,Recent3Y,Recent5Y,Y2018,Y2017,Y2016,Y2015,Y2014))

fo.close()
