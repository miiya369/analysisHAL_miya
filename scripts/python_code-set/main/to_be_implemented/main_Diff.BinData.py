#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from MiscFuncs.IO_BinData     import input_bin_data

iFname1 = None
iFname2 = None

###### Functions for arguments
def usage(ARGV):
    print
    print "usage: %s [ifname 1] [ifname 2]" % ARGV[0]
    print; quit()

def set_args(ARGC, ARGV):
    global iFname1, iFname2
    
    iFname1 = ARGV[1].strip()
    iFname2 = ARGV[2].strip()
    
###### Main part
if __name__ == "__main__":
    from numpy import empty
    
    argv = sys.argv; argc = len(argv)
    
    if (argc != 3):
        usage(argv)
    
    set_args(argc, argv)
    
    yData1, xData1, Dummy = input_bin_data(iFname1)
    if (xData1 is yData1 is Dummy is None):
        quit()
    yData2, xData2, Dummy = input_bin_data(iFname2)
    if (xData2 is yData2 is Dummy is None):
        quit()
    
    if (len(xData1) != len(xData2)):
        print("ERROR: The lengths of each file are differ, exit.\n"); quit()
    
    Ndata = len(xData1)
    
    for i in range(Ndata):
        if (xData1[i] != xData2[i]):
            print("ERROR: The x-data of each file are differ, exit.\n"); quit()
    
    Nconf = len(yData1[0, :])
    tmp_d = empty(Nconf)
    
    for iconf in range(Nconf):
        tmp_d[iconf] = 0.0
    
    for idata in range(Ndata):
        for iconf in range(Nconf):
            tmp_d[iconf] =     yData1[idata, iconf] - yData2[idata, iconf]
            #tmp_d[iconf] += abs(yData1[idata, iconf] - yData2[idata, iconf])
        mean, err = make_mean_err(tmp_d)
        print("%lf %1.16lf %1.16lf" % (xData1[idata], mean, err))
    
    mean, err = make_mean_err(tmp_d)
    print("#Dif: %1.16lf +/- %1.16lf" % (mean, err))
