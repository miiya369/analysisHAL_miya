#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from MiscFuncs.IO_BinData     import input_bin_data

iFname = None

###### Functions for arguments
def usage(ARGV):
    print
    print "usage: %s [bin-data]" % ARGV[0]
    print
    quit()

def set_args(ARGC, ARGV):
    global iFname
    
    iFname = ARGV[1].strip()

###### Main part
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc != 2):
        usage(argv)
    
    set_args(argc, argv)
        
    yData, xData, eData = input_bin_data(iFname)
    if (xData is yData is eData is None):
        quit()
    
    for idata in range(len(xData)):
        mean, err = make_mean_err(yData[idata, :])
        print("%lf %1.16e %1.16e" % (xData[idata], mean, err)),
        for iconf in range(len(yData[idata, :])):
            print(" %1.16e" % yData[idata, iconf]),
        print
