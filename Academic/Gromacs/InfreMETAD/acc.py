#!/usr/bin/python
from math import *
from itertools import ifilterfalse


kbt=2.496 
Temperature=300   # Temperature in Kelvins
dtmd=1000           # timestep in COLVAR file (fs)
metabiascol=4      # Column of metadynamics bias (index starting from 1)
fname='inmetad'     # name of COLVAR file


metabiascol=metabiascol-1  # python index starts from 0

conversion=(1e-12)

def iscomment(s):
    return s.startswith('#')



for i in [2,3,4,5,6,7,8,9,10,11,13,15,16,17,20]:         #
   #print i+1
   #print '----------------------------------'
   string=fname+'_'+str(i) + '/' + 'COLVAR_NEW'
   accfile=open('inmetad_'+str(i)+'/'+'acc.dat', 'w')
   f=open(string,'r')

   timesum=0.0
   countlines=0.0
   for line in ifilterfalse(iscomment,f):
      line=line.strip()
      columns=line.split()
      bias=float(columns[metabiascol])
      timesum=timesum+exp(300/kbt/Temperature*bias)*dtmd
      countlines= countlines+1
      accfactor = timesum/(dtmd*countlines)
      if countlines%10 == 0:
         accfile.write(str(countlines)+' '+str(accfactor)+' ' + str(countlines*accfactor)+'\n')
     # if float(columns[water])>=waterEnd and float(columns[dis])>=disEnd:
     # 	print timesum*conversion
     # 	break;
