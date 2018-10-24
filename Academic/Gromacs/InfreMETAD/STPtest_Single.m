% module add matlab/r2016a
% matlab < STPtest_Single.m >STPtest_Single.log

load('times.txt');

H=STPtest(times,min(times)./10,max(times).*10,1E4,0.05)
