# -*- coding: utf-8 -*-

import math
KbT = 2.496           # KbT at room temperature
dt = 0.2              # Time in picosecond (ps) between two lines of COLVAR
sum_bias = 0

############Remove the overlapping lines of COLVAR##########
i = 0
time0 = -1
fo = open('COLVAR-no-overlap', 'w')
with open('COLVAR', 'r') as fi:
    for line in fi:
        # Only write the first comment lines, omit the commet lines for restarting runs.
        if line.startswith('#') and i < 3:
            fo.write(line)
            i += 1
        elif not line.startswith('#'):
            time = int(line.strip().split()[0].split('.')[0])
            # Coninuous output
            if time > time0:
                time0 = time
                fo.write(line)
            # Redundant output, ignore it.
            elif time < time0:
                continue
        
fo.close()

###############Process the non-overlapping COLVAR##########
with open('COLVAR-no-overlap', 'r') as fi:
    for line in fi:
        if not line.startswith('#'):
            terms = line.split()
            time = float(terms[0])      # The instant metadynamics simulation time (ps).
            dist = float(terms[1])
            bias = float(terms[3])
            if dist < 2.0:
                sum_bias += math.exp(1 / KbT * bias)
            else:
                # When the distance is bigger than the cut-off, it is seen as dissociated.
                metad_time = time       # The metadynamics simulation (ps) when the ligand dissociates.
                break
    
    resid_time = dt * sum_bias          # Residence time in ps.
    alpha = resid_time / metad_time     # Both in ps.
    meta_time = metad_time * math.pow(10, -3)          # Convet ps to ns
    resid_time = dt * sum_bias * math.pow(10, -9)      # Convert ps to ms.

print "Metadynamics time: {:8.3f} ns".format(metad_time)
print "Residence time: {:8.3f} ms".format(resid_time)
print "Alpha is: {:8.3f}".format(alpha)