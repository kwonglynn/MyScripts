# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 19:22:37 2019

@author: Guanglin Kuang
"""

import tushare as ts
import datetime

def cal_com_interest(fund_changes):
    # This is the right way to calculate fund increae.
    list_changes = list(fund_changes)
    com_interest = 1.0
    for change in list_changes:
        # For dividend
        if abs(change) > 11:
            print (change)
            change = 0
        com_interest = com_interest * (1.0 + change/100.0)
    increase = (com_interest - 1.0) * 100
    return increase

# 循环处理每个基金
def process_fund(code, fo):
    # By defalut only the values in one year are returned.
    fund = ts.fund.nav.get_nav_history(code, start=start_date, end=end_date)
    changes = fund.change.dropna()

    #Return in one, two and three years
    this_year = cal_com_interest(changes.loc[today:'2019-01-01'])
#    ThisYear = changes.loc[today:'2019-01-01'].sum()
    
    one_year_ago = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365)).strftime('%Y-%m-%d')
    recent1Y = cal_com_interest(changes.loc[today:one_year_ago])
#    Recent1Y = changes.loc[today:OneYearAgo].sum()

    two_years_ago = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 2)).strftime('%Y-%m-%d')
#    Recent2Y = changes.loc[today:TwoYearsAgo].sum()
    recent2Y = cal_com_interest(changes.loc[today:two_years_ago])

    three_years_ago = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 3)).strftime('%Y-%m-%d')
#    Recent3Y = changes.loc[today:ThreeYearsAgo].sum()
    recent3Y = cal_com_interest(changes.loc[today:three_years_ago])

    five_years_ago = (datetime.datetime.strptime(today, '%Y-%m-%d') - datetime.timedelta(days = 365 * 5)).strftime('%Y-%m-%d')
#    Recent5Y = changes.loc[today:FiveYearsAgo].sum()
    recent5Y = cal_com_interest(changes.loc[today:five_years_ago])
    
    Y2018 = cal_com_interest(changes.loc['2018-12-31':'2018-01-01'])
    Y2017 = cal_com_interest(changes.loc['2017-12-31':'2017-01-01'])
    Y2016 = cal_com_interest(changes.loc['2016-12-31':'2016-01-01'])
    Y2015 = cal_com_interest(changes.loc['2015-12-31':'2015-01-01'])
    Y2014 = cal_com_interest(changes.loc['2014-12-31':'2014-01-01'])     

    fo.write("'%s'\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n" % \
             (code, this_year, recent1Y, recent2Y, recent3Y, recent5Y, Y2018, Y2017, Y2016, Y2015, Y2014))

##################################Main#######################################
today = datetime.date.today().strftime('%Y-%m-%d')

# The start date should be 6 years ago!
start_date = '2013-01-01'
end_date = '2019-03-21'

fund_types = ['flexible', 'large', 'small']

# General funds
# 1. Flexible funds
if 'flexible' in fund_types:
    fo = open('Flexible-Fund_Rating-%s.txt' % today, 'w')
    fo.write("Code\tThisYear\tRecent1Y\tRecent2Y\tRecent3Y\tRecent5Y\tY2018\tY2017\tY2016\tY2015\tY2014\n")
    fund_list = ['000742', '001102', '020005', '450002', '000628', 
                 '001208', '000577', '519697', '519688', '000527', '519091']
    for code in fund_list:
        process_fund(code, fo)
    fo.close()


# 2. Large funds
if 'large' in fund_types:
    fo = open('Large-Fund_Rating-%s.txt' % today, 'w')
    fo.write("Code\tThisYear\tRecent1Y\tRecent2Y\tRecent3Y\tRecent5Y\tY2018\tY2017\tY2016\tY2015\tY2014\n")
    fund_list = ['040008', '162605', '260108', '000619', '110011', 
                 '180012', '540012', '540006', '519068', '519066']
    for code in fund_list:
        process_fund(code, fo)
    fo.close()

# 3. Small funds
if 'small' in fund_types:
    fo = open('Small-Fund_Rating-%s.txt' % today, 'w')
    fo.write("Code\tThisYear\tRecent1Y\tRecent2Y\tRecent3Y\tRecent5Y\tY2018\tY2017\tY2016\tY2015\tY2014\n")
    fund_list = ['000751', '020026', '090013', '001938', '519732',
                 '001985', '001178', '519712', '519736', '001410',
                 '040011', '420003']
    for code in fund_list:
        process_fund(code, fo)
    fo.close()

