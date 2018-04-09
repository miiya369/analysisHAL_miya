#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from MiscFuncs.IO_BinData     import input_bin_data

iFname = None
r_max  = 100.0

###### Functions for arguments
def usage(ARGV):
    print
    print "usage: %s [ifname] {[r max]}" % ARGV[0]
    print; quit()

def set_args(ARGC, ARGV):
    global iFname, r_max
    
    iFname = ARGV[1].strip()
    if (ARGC > 2):
        r_max  = float(ARGV[2].strip())
    
###### Main part
if __name__ == "__main__":
    from numpy import empty
    
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1 or argc > 3):
        usage(argv)
    
    set_args(argc, argv)
    print("# r_max = %lf" % r_max)
    
    yData, xData, Dummy = input_bin_data(iFname)
    if (xData is yData is Dummy is None):
        quit()
    
    Nconf = len(yData[0, :])
    Ndata = len(yData[:, 0])
    
    r2V = empty(Nconf)
    
    for iconf in range(Nconf):
        r2V[iconf] = 0.0
        for idata in range(Ndata):
            if (r_max > xData[idata]):
                r2V[iconf] += yData[idata, iconf] * xData[idata]**2
    
    mean, err = make_mean_err(r2V)
    print("r2V: %1.16lf +/- %1.16lf" % (mean, err))
