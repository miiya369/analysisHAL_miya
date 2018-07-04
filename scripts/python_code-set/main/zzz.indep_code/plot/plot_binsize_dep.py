#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc, patches

read_colm = 1
tmin      = 9
tmax      = 12

ofname = None
#ofname = "/Users/miiya/Desktop/effmass_kaon.pdf"

### =========================== Main =========================== ###

def main(ifpaths):
    Corr = read_correlators(ifpaths[0])
    for i in range(1, len(ifpaths)):
        Corr = np.append(Corr, read_correlators(ifpaths[i]), axis=0)
    Nfile = len(Corr[:,0,0])
    Ntime = len(Corr[0,:,0])
    Ncolm = len(Corr[0,0,:])
    print("# Read correlators: Nfile=%d, Ntime=%d, Ncolm=%d" % (Nfile, Ntime, Ncolm))
    
    bsize  = []
    jkerr  = []
    jkerer = []
    for ibsize in range(1, Nfile):
        if (Nfile % ibsize != 0):
            continue
        jkCorr = make_jk_samples(Corr[:,tmin:tmax+2,read_colm], ibsize)
        jkMass = np.array([np.log(jkCorr[:,it]/jkCorr[:,it+1]) for it in range(tmax-tmin+1)])
        jkErrs = np.std(jkMass, axis=1) * np.sqrt(len(jkMass[0,:])-1)
        
        bsize .append(ibsize)
        jkerr .append(np.mean(jkErrs))
        jkerer.append(np.std (jkErrs))
    
    # --- plot --- #
    rc('text' , usetex=True)
    rc('font' ,**{'family'    : 'Times New Roman',
                  'weight'    : 'bold',
                  'size'      : 11})
    rc('xtick',**{'labelsize' : 12})
    rc('ytick',**{'labelsize' : 12})
    rc('axes' ,**{'labelsize' : 12})
    rc('axes' ,**{'linewidth' : 2})
    
    plt.ion(); fig = plt.figure()
    fig.patch.set_facecolor('white')
    axe = fig.add_axes((0.15, 0.15, 0.8, 0.8))
    axe.patch.set_facecolor('white')
    
    axe.errorbar(np.array(bsize), np.array(jkerr), yerr=(np.array(jkerer)), fmt='none', #fmt=None,
                 marker=None, lw=3, ecolor='red', capthick=3, capsize=3)
    
    axe.grid(which='major',color='gray',linestyle='--')
    axe.set_ylabel('Jack-knife error for effective mass')
    axe.set_xlabel('bin-size')
    axe.yaxis.set_label_coords(-0.12,0.5)
    axe.xaxis.set_label_coords(0.5,-0.1)
    
    if (ofname is not None):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    else:
        plt.show(); print("# push Enter to end"); input('# ')

    return 0

### ============================================================ ###

def read_correlators(ifpath, auto_fb = True):
    from glob    import glob, iglob
    from os.path import dirname, basename
    Corr        = np.array([np.loadtxt(ifname) for ifname in iglob(ifpath)])
    ifpath_anti = dirname(ifpath) + "/anti" + basename(ifpath)
    if (len(glob(ifpath_anti)) == 0 or not auto_fb):
        return Corr
    else:
        print("# Read anti-hadron...")
        Corr_anti = np.array([np.loadtxt(ifname) for ifname in iglob(ifpath_anti)])
        return (Corr + np.roll(Corr_anti[:,::-1,:], 1, axis=1)) / 2.0

def make_jk_samples (idata, Bsize):
    Ndata = len(idata[:,0])
    Nbin  = Ndata // Bsize
    return np.array([(np.sum(idata, axis=0)-np.sum(idata[i*Bsize:(i+1)*Bsize, :], axis=0))
                     for i in range(Nbin)]) / float(Ndata-Bsize)

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc < 2):
        exit("usage: python %s [ifile(s)1] [ifile(s)2] ..." % basename(argv[0]))
    
    if (main(argv[1:]) != 0):
        exit("ERROR EXIT.")
