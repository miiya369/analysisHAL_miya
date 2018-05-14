#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from MiscFuncs.IO_BinData     import input_bin_data, output_bin_data

Nadd   = 0
iFname = []
oFname = None
factor = []

###### Functions for arguments
def usage(ARGV):
    print
    print "usage: %s [ofname] [factor1] [ifname1] [factor2] [ifname2] ..." % ARGV[0]
    print
    print " Note: The factor is not supposed to contain any operators (except to the sign)."
    print
    quit()

def set_args(ARGC, ARGV):
    global iFname, oFname, Nadd, factor
    
    oFname = ARGV[1].strip()
    Nadd   = (ARGC - 2) / 2
    
    for iadd in range(2, ARGC, 2):
        factor.append(float(ARGV[iadd+0])       )
        iFname.append(      ARGV[iadd+1].strip())
    
###### Main part
if __name__ == "__main__":
    from numpy import empty
    
    argv = sys.argv; argc = len(argv)
    
    if (argc < 6 or argc % 2 != 0):
        usage(argv)
    
    set_args(argc, argv)
    
    yData, xData, Dummy = input_bin_data(iFname[0])
    if (xData is yData is Dummy is None):
        quit()
    
    Nconf = len(yData[0, :])
    Ndata = len(yData[:, 0])
    
    yData *= factor[0]
    
    for iadd in range(1, Nadd):
        tmpyData, tmpxData, Dummy = input_bin_data(iFname[iadd])
        if (tmpxData is tmpyData is Dummy is None):
            quit()
        
        for idata in range(Ndata):
            if (xData[idata] != tmpxData[idata]):
                print("\nERROR: xData is differ, exit.\n"); quit()
            for iconf in range(Nconf):
                yData[idata, iconf] += (factor[iadd] * tmpyData[idata, iconf])
        del tmpyData, tmpxData
    
    output_bin_data(oFname, yData, xData)
