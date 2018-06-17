#!/afs/pdc.kth.se/home/g/guanglin/nobackup/software/anaconda2/bin/python
##Change the path to python for your own system!!!

"""
Author:
    Guanglin Kuang <guanglin@kth.se>

Usage:
    calculate_rmsd_dist.py [options]

Options:
    -r, --ref <ref>             The reference structure [default: complex_box.gro].
    -t, --top <top>      	The topology file for the trajectory [default: complex_box.gro].
    -j, --traj <traj>           The trajectory file after PBC treatment [default: md_noPBC.xtc].
    -n, --N_frames <N_frames>  	The number frames to process.
    -p, --prefix <prefix>	The prefix for the output files [default: MDAnalysis].

Note:
    1, The toplogy file for Gromacs format is the gro file, not the top file.

Example:
    ./calculate_rmsd_dist.py -r complex_box.gro -t complex_box.gro -j scaled_noPBC.xtc -p MDAnalysis

"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
from docopt import docopt

def trans_rot(ref, mobile):
    mobile0 = mobile.select_atoms('backbone').positions - mobile.atoms.center_of_mass()
    ref0 = ref.select_atoms('backbone').positions - ref.atoms.center_of_mass()

    R, rmsd = align.rotation_matrix(mobile0, ref0)

    mobile.atoms.translate(-mobile.select_atoms('backbone').center_of_mass())
    mobile.atoms.translate(ref.select_atoms('backbone').center_of_mass())
    mobile.atoms.rotate(R)

def measure_rmsd(ref, mobile):
    backbone_ref = ref.select_atoms("backbone")
    LIG_ref = ref.select_atoms("resname LIG")

    backbone_mobile = mobile.select_atoms("backbone")
    LIG_mobile = mobile.select_atoms("resname LIG")

    rmsd_backbone = rms.rmsd(backbone_mobile.positions, backbone_ref.positions)
    rmsd_LIG = rms.rmsd(LIG_mobile.positions, LIG_ref.positions)

    return rmsd_backbone, rmsd_LIG

def measure_dist(ref, mobile):
    COM_LIG = mobile.select_atoms("resname LIG").center_of_mass()
    COM_poc = ref.select_atoms("protein and around 6 resname LIG").center_of_mass()
    dist = np.linalg.norm(COM_LIG - COM_poc)
    return dist

##########Input Stop here###############################
if __name__ == '__main__':
    ####Treat the options####
    opts = docopt(__doc__)
    print opts
    ref = mda.Universe(opts['--ref'])
    mobile = mda.Universe(opts['--top'], opts['--traj'])

    if opts['--N_frames']:
        N_frames = int(opts['--N_frames'])
    else:
	N_frames = mobile.trajectory.n_frames

    prefix = opts['--prefix']

    data = []
    for ts in mobile.trajectory[0:N_frames]:
        time = ts.time/1000
        trans_rot(ref, mobile)
        rmsd_backbone, rmsd_LIG = measure_rmsd(ref, mobile)
        dist = measure_dist(ref, mobile)
        data.append([time, rmsd_backbone, rmsd_LIG, dist])

    fo = open('%s.dat' % prefix, 'w')
    fo.write("#Time\tRMSD_backbone\tRMSD_LIG\tCOM_distance\n")

    for item in data:
#        fo.write('\t'.join(format("%8.2f" % i) for i in item) + '\n')
        fo.write('{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n'.format(*item))

    fo.close()

    data = np.array(data).transpose()

    ax1 = plt.subplot(111)
    ax1.plot(data[0], data[1], label='RMSD backbone')
    ax1.plot(data[0], data[2], label='RMSD ligand')
    ax1.plot(data[0], data[3], label='COM distance')
    ax1.set_xlabel("Time (ps)")
    ax1.set_ylabel(r'$\AA$')
    ax1.legend(loc='best')
    plt.savefig('%s.png' % prefix, dpi=600)

