#!/usr/bin/python

from __future__ import print_function

import numpy as np

Tsize = 32

def main (wdir):
    from glob import glob, iglob
    
    for iconf in iglob(wdir+'/*'):
        conf_base = basename(iconf)
        print("Averaging in "+conf_base+" ...")
        for ihad in iglob(iconf+"/*.+000.000.000.000."+conf_base):
            had_name = basename(ihad).split('.')[0]
            ifcorrs  = glob(iconf+"/"+had_name+".+0??.000.000.000."+conf_base)
            ofcorr   =      iconf+"/"+had_name+".+Ave.000.000.000."+conf_base
            if (len(ifcorrs) != Tsize):
                print("ERROR in %s: Unexpected #.ave." % (iconf+"/"+had_name)); return -1            
            average_textdata(ofcorr, ifcorrs)
            print("   DONE (#.ave=%d): " % len(ifcorrs) +had_name)
    return 0

def average_textdata (ofpath, ifpaths):
    ave_data = np.mean(np.array([np.loadtxt(ifpaths[i]) for i in range(len(ifpaths))]), axis=0)
    np.savetxt(ofpath, ave_data, fmt='%4d'+' %1.16e'*(len(ave_data[0, :])-1))
    return 0

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename

    if (argc != 2):
        exit("usage: python %s [wdir]" % basename(argv[0]))
    
    if (main(argv[1]) != 0):
        exit("ERROR EXIT.")
