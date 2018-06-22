#!/afs/pdc.kth.se/home/g/guanglin/nobackup/software/anaconda2/bin/python
##Change the path to python for your own system!!!

"""
Author:
    Guanglin Kuang <guanglin@kth.se>

Usage:
    Extract_Umbrella_Frames.py [options]

Options:
    -r, --ref <ref>             The reference structure [default: complex_box.gro].
    -t, --top <top>             The topology file for the trajectory [default: complex_box.gro].
    -j, --traj <traj>           The trajectory file after PBC treatment [default: pull_noPBC.xtc].
    -n, --N_frames <N_frames>   The number frames to process, default is to process frames.
    -p, --prefix <prefix>       The prefix for the output files [default: Umbrella_Distances.dat].

Note:
    1, The toplogy file for Gromacs format is the gro file, not the top file.
    2, The PBC of the trajectory should be processed:
    gmx_seq trjconv -s pull.tpr -f pull.xtc -pbc mol -ur compact -n index.ndx -o pull_noPBC.xtc
    
Example:
    python Extract_Umbrella_Frames.py -r complex_box.gro -t complex_box.gro -j pull_noPBC.xtc -p Umbrella_Distances

"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
from docopt import docopt
import math

def trans_rot(ref, mobile):
    mobile0 = mobile.select_atoms('backbone').positions - mobile.atoms.center_of_mass()
    ref0 = ref.select_atoms('backbone').positions - ref.atoms.center_of_mass()

    R, rmsd = align.rotation_matrix(mobile0, ref0)

    mobile.atoms.translate(-mobile.select_atoms('backbone').center_of_mass())
    mobile.atoms.translate(ref.select_atoms('backbone').center_of_mass())
    mobile.atoms.rotate(R)

def measure_dist(ref, mobile):
    COM_LIG = mobile.select_atoms("resname LIG").center_of_mass()
    COM_poc = ref.select_atoms("backbone and around 6 resname LIG").center_of_mass()
    dist = np.linalg.norm(COM_LIG - COM_poc)
    return dist

##########Input Stop here###############################
if __name__ == '__main__':
    ####Treat the options####
    opts = docopt(__doc__)
    print (opts)
    ref = mda.Universe(opts['--ref'])
    mobile = mda.Universe(opts['--top'], opts['--traj'])

    if opts['--N_frames']:
        N_frames = int(opts['--N_frames'])
    else:
        N_frames = mobile.trajectory.n_frames

    prefix = opts['--prefix']

    data = []
    dists_want = range(1,21)
    dists_now = []
    prot_LIG = mobile.select_atoms("protein or resname LIG")
    for ts in mobile.trajectory[0:N_frames]:
        if ts.frame % 1000 == 0:
            print ("Processing frame: {:d}".format(ts.frame))
        time = ts.time/1000
        trans_rot(ref, mobile)
        dist = measure_dist(ref, mobile)
        dist_int = int(math.floor(dist))
        if (dist_int in dists_want) and (dist_int not in dists_now):
            dists_now.append(dist_int)
            data.append([time, dist, 'True'])
            prot_LIG.write("complex_{:d}.gro".format(dist_int))          
        else:
            data.append([time, dist, ''])

    fo = open('%s.dat' % prefix, 'w')
    fo.write("#Time\tCOM_distance\tSelect\n")

    for item in data:
        fo.write('{:.2f}\t{:.2f}\t{:s}\n'.format(*item))

    fo.close()

