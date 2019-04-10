import pandas as pd
import numpy as np

## Get the number of lines and then determine the number of files:
df = pd.read_csv('HILLS', comment='#')
N = 1 + len(df.index) // 100000

## Generate a list of file objects to write the separated HILLS.
fos = []
for i in range(N):
    file_name = 'HILLS-' + str(i)
    fo = open(file_name, 'w')
    fos.append(fo)

fi = open('HILLS', 'r')
for line_no, line in enumerate(fi):
    if line.startswith('#'):
        if line_no < 10:
            for fo in fos:
                fo.write(line)
            continue
        else:        
            continue
    time = int(line.split()[0].split('.')[0])
    n = time // 100000
    # There is overlap in each HILLS file.
    for i in range(n, N):
        fos[i].write(line)

fi.close()
for fo in fos:
    fo.close()

## 1. Check the line numbers in each HILLS:
## find . -name  "HILLS*" | sort -V | xargs -i wc -l {}

## 2. Run plumed on each HILLS file
## for i in {0..56}; do plumed sum_hills --hills HILLS-$i --bin 98,99 --outfile fes-$i.dat; done

## 3. Run FES_METAD.py
## for i in {0..56}; do python FES_METAD.py fes-$i.dat -o fes-$i.png

## 4. Calculate the similarity between the FES images.
## python compare_FES_SSIM.py
