fos = []
for i in range(1, 7):
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
    N = time // 1000000
    # There is overlap in each HILLS file.
    for i in range(N, 6):
        fos[i].write(line)

fi.close()
for fo in fos:
    fo.close()

## Plumed commands:
## for i in {1..6};do plumed sum_hills --hills HILLS-$i --bin 98,99 --outfile fes-$i.dat; done
