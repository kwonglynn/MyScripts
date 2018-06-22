# -*- coding: utf-8 -*-
fo = open('Scaled_times.dat' , 'w')

with open('MDAnalysis_path.dat', 'r') as fi_path:
    for path in fi_path:
        path = path.strip()
        with open(path, 'r') as fi_data:
            for data in fi_data:
                terms = data.split()
                time = float(terms[0])
                dist = float(terms[3])
                if dist >= 2.0:
                    fo.write("{.2f}\n".format(time))
                    
fo.close()
                
    