#!/bin/sh
#

file=Raw_times.txt
awk '!/#/ {print $5}' $file >times.txt
#/home/zy/MATLAB/R2016a/bin/matlab -nojvm < error_Nsim.m >SD_Nsim.log
/home/zy/MATLAB/R2016a/bin/matlab < error_Nsim.m >SD_Nsim.log
cat SD_Nsim.txt
