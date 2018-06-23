# !/usr/bin/python
# -*- coding: utf-8 -*-
import datetime
import tushare as ts
import pandas as pd
import numpy as np


def CalcDownside(returns):
    #mu = 0.35 / 244 * 100
    mu = returns.mean()
    downs = returns[returns < mu]
    downside = 100 * (sum((mu-downs)**2) / len(returns)) ** 0.5
    downside = float("%4.1f" % downside)
    return downside

def CalcYear(code, StartDate, EndDate):
    # 2016-2017
    fund1 = ts.fund.nav.get_nav_history(code, start='2016-02-01', end='2017-02-01')
    price1 = fund1.total.dropna()
    year1 = (price1[0] - price1[-1]) / price1[-1] * 100
    year1 = int(year1)
    
    # 2017-2018
    fund2 = ts.fund.nav.get_nav_history(code, start='2017-01-01', end='2018-01-01')
    price2 = fund2.total.dropna()
    year2 = (price2[0] - price2[-1]) / price2[-1] * 100
    year2 = int(year2)
    
    # The past one year from now.
    tEndDate = datetime.datetime.now()
    tStartDate = tEndDate + datetime.timedelta(days=-365)
    sEndDate = datetime.datetime.strftime(tEndDate, '%Y-%m-%d')
    sStartDate = datetime.datetime.strftime(tStartDate, '%Y-%m-%d')
    fund = ts.fund.nav.get_nav_history(code, start=sStartDate, end=sEndDate)
    price = fund.total.dropna()
    year = (price[0] - price[-1]) / price[-1] * 100
    year = int(year)

    years = [year1, year2, year]
    
    if year > 35 and np.mean(years) > 30:
        keep = 'True'
    else:
        keep = 'False'
        
    return keep, year, years

def CalcSeason(code, StartDate, EndDate):
    tStartDate = datetime.datetime.strptime(StartDate, '%Y-%m-%d')
    tEndDate = datetime.datetime.strptime(EndDate, '%Y-%m-%d')
    
    seasons = []
    date1 = StartDate                   # date1 is string type
    date2 = tStartDate + datetime.timedelta(days=90)    # date2 is datetime type.
    while date2 <= tEndDate:
        # Process the first month.
        date2 = datetime.datetime.strftime(date2, '%Y-%m-%d')   # datetime to string    
        fund = ts.fund.nav.get_nav_history(code, start=date1, end=date2)
        prices = fund.total.dropna()
        increase =  (prices[0]-min(prices)) /min(prices) * 100
        increase = int(increase)
        seasons.append(increase)
            
        # Prepare the next month
        date1 = date2
        date2 = datetime.datetime.strptime(date2, '%Y-%m-%d')   # string to datetime
        date2 = date2 + datetime.timedelta(days=90)
    return seasons

def CalcDrawdown(code, StartDate, EndDate):
    tStartDate = datetime.datetime.strptime(StartDate, '%Y-%m-%d')
    tEndDate = datetime.datetime.strptime(EndDate, '%Y-%m-%d')
    
    MaxDrawdowns = []
    date1 = StartDate                   # date1 is string type
    date2 = tStartDate + datetime.timedelta(days=30)    # date2 is datetime type.
    while date2 <= tEndDate:
        # Process the first month.
        date2 = datetime.datetime.strftime(date2, '%Y-%m-%d')   # datetime to string    
        fund = ts.fund.nav.get_nav_history(code, start=date1, end=date2)
        prices = fund.total.dropna()
        low = min(prices[0:5])
        high = max(prices)
        drawdown = (low - high) / high * 100
        drawdown = int(drawdown)
        MaxDrawdowns.append(drawdown)
        
        # Prepare the next month
        date1 = date2
        date2 = datetime.datetime.strptime(date2, '%Y-%m-%d')   # string to datetime
        date2 = date2 + datetime.timedelta(days=30)
    
    AveDrawdown = np.mean(MaxDrawdowns)
    MaxDrawdown = min(MaxDrawdowns)
    return AveDrawdown, MaxDrawdown

##################################Main########################################
StartDate = '2016-02-01'
EndDate = '2018-01-20'
 
test = 'True'
if test == 'True': 
    codes = ['110022','540006','590008','519606','090013','020026']

fo = open('Fund_Mining2.dat','w')
fo.write("Code\tYear%\tYears%\tSTD\tDownside%\tAveDrawdown%\tMaxDrawdown%\tSeasons%\n")

if test != 'True':
    with open('EquityFunds.dat','r') as fi:
        codes=fi.readlines()
    
for code in codes:
    code = code.strip()
    fund = ts.fund.nav.get_nav_history(code, start=StartDate, end=EndDate)
    prices = fund.total.dropna()
    returns = fund.change.dropna()
    
    if len(prices) < 480:
        continue
        
    keep, year, years = CalcYear(code, StartDate, EndDate)
    if keep != 'True': 
        continue
    
    print code    
    std = returns.std()
    downside = CalcDownside(returns)
    AveDrawdown, MaxDrawdown = CalcDrawdown(code, StartDate, EndDate)
    seasons =  CalcSeason(code, StartDate, EndDate)

    fo.write("'%s'\t%d\t%s\t%4.2f\t%4.1f\t%4.1f\t%4.1f\t%s\n" % \
             (code, year, str(years), std, downside, AveDrawdown, MaxDrawdown, str(seasons)))
    
fo.close()
        

    
    
    

