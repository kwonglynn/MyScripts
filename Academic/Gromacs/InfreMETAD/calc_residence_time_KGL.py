# -*- coding: utf-8 -*-
'''
Created on Sept. 13, 2018

@author: Guanglin Kuang (guanglin@kth.se)
'''

import math

############Remove the overlapping lines of COLVAR##########
def ProcessCOLVAR(COLVAR_in, COLVAR_out):
    i = 0
    time0 = -1
    fo = open(COLVAR_out, 'w')
    with open(COLVAR_in, 'r') as fi:
        for line in fi:
            # Only write the first comment lines, omit the commet lines for restarting runs.
            if line.startswith('#') and i < 3:
                fo.write(line)
                i += 1
            elif not line.startswith('#') and len(line) > 1:
                time = int(line.strip().split()[0].split('.')[0]) # Only lines with interger ps times are kept.
            # Coninuous output
                if time > time0:
                    time0 = time
                    fo.write(line)
            # Redundant output, ignore it.
                else:
                    continue
        
    fo.close()

###############Process the non-overlapping COLVAR##########
def CalMetad(COLVAR, KbT, dt):
    sum_bias = 0
    with open(COLVAR, 'r') as fi:
        for line in fi:
            if not line.startswith('#'):
                terms = line.split()
                time = float(terms[0])      # The instant metadynamics simulation time (ps).
                dist = float(terms[1])
                bias = float(terms[3])
                if dist < 1.9:
                    sum_bias += math.exp(1 / KbT * bias)
                else:
                    # When the distance is bigger than the cut-off, it is seen as dissociated.
                    metad_time = time       # The metadynamics simulation (ps) when the ligand dissociates.
                    break
    
        resid_time = dt * sum_bias          # Residence time in ps.
        alpha = resid_time / metad_time     # Both in ps.
        metad_time = metad_time * math.pow(10, -3)          # Convet ps to ns
        resid_time = resid_time * math.pow(10, -9)      # Convert ps to ms.
        
    return metad_time, resid_time, alpha

def main():
    ###Define the number of simulations (N), the temperature (T), and time step of COLVAR (dt). Note dt is NOT the time step of MD simulation.
    N = 20                  # Number of independent simulations.
    T = 300                 # The temperature of the simulation in K
    dt = 1.0                # Time in picosecond (ps) between two lines of COLVAR after removing overlapping lines.
    ###
    KbT0 = 2.479            # KbT at 298 K, in kJ/mol
    KbT = T / 298 * KbT0    # Calculate the KbT at the simulation temperature

    fo = open("InfreMetad_Results.dat", 'w')
    fo.write("#Metad_time(ns)\tResid_time(ms)\tAlpha")
    
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20]:
    	#for i in range(1, N+1):
	print i
        ProcessCOLVAR("Run%d/COLVAR" % i, "Run%d/COLVAR-no-overlap" % i)
        metad_time, resid_time, alpha = CalMetad("Run%d/COLVAR-no-overlap" % i, KbT, dt)
        fo.write("{:8.3f}\t{:8.2f}\t{:.3e}\n".format(metad_time, resid_time, alpha))
    
    fo.close()
    
if __name__ == "__main__":
    main()
