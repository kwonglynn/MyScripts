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

# The start date should be 6 years ago!
StartDate = '2013-01-01'
EndDate = '2019-03-03'

AllFundsList = ['519068', '162605', '000628', '000527', '160505', '450002', '020026', '090013', '519069', '519066', \
                '000577', '110011', '180012', '040008', '519732', '000619', '519736', '260108', '540012', '001938', \
                '000742', '000751', '340008', '540006', '070032']

def CalComInterest(InputChanges):
    # This is the right way to calculate fund increae.
    ListChanges = list(InputChanges)
    ComInterest = 1.0
    for change in ListChanges:
        # For dividend
        if abs(change) > 10:
            print (change)
            change = 0
        ComInterest = ComInterest * (1.0 + change/100.0)
    increase = ComInterest - 1.0
    return increase

# 循环处理每个基金
for code in AllFundsList:
    print (code)
    # By defalut only the values in one year are returned.
    fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
    changes = fund.change.dropna()

    #Return in one, two and three years
    ThisYear = CalComInterest(changes.loc[today:'2019-01-01'])
#    ThisYear = changes.loc[today:'2019-01-01'].sum()
    
    OneYearAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365)).strftime('%Y-%m-%d')
    Recent1Y = CalComInterest(changes.loc[today:OneYearAgo])
#    Recent1Y = changes.loc[today:OneYearAgo].sum()

    TwoYearsAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 2)).strftime('%Y-%m-%d')
#    Recent2Y = changes.loc[today:TwoYearsAgo].sum()
    Recent2Y = CalComInterest(changes.loc[today:TwoYearsAgo])

    ThreeYearsAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 3)).strftime('%Y-%m-%d')
#    Recent3Y = changes.loc[today:ThreeYearsAgo].sum()
    Recent3Y = CalComInterest(changes.loc[today:ThreeYearsAgo])

    FiveYearsAgo = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 5)).strftime('%Y-%m-%d')
#    Recent5Y = changes.loc[today:FiveYearsAgo].sum()
    Recent5Y = CalComInterest(changes.loc[today:FiveYearsAgo])

#    Y2018 = changes.loc['2018-12-31':'2018-01-01'].sum()
#    Y2017 = changes.loc['2017-12-31':'2017-01-01'].sum()
#    Y2016 = changes.loc['2016-12-31':'2016-01-01'].sum()
#    Y2015 = changes.loc['2015-12-31':'2015-01-01'].sum()
#    Y2014 = changes.loc['2014-12-31':'2014-01-01'].sum()    
    
    Y2018 = CalComInterest(changes.loc['2018-12-31':'2018-01-01'])
    Y2017 = CalComInterest(changes.loc['2017-12-31':'2017-01-01'])
    Y2016 = CalComInterest(changes.loc['2016-12-31':'2016-01-01'])
    Y2015 = CalComInterest(changes.loc['2015-12-31':'2015-01-01'])
    Y2014 = CalComInterest(changes.loc['2014-12-31':'2014-01-01'])     

    fo.write("'%s'\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n" % \
             (code,ThisYear,Recent1Y,Recent2Y,Recent3Y,Recent5Y,Y2018,Y2017,Y2016,Y2015,Y2014))

fo.close()
