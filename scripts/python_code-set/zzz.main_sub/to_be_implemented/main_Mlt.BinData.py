#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from MiscFuncs.IO_BinData     import input_bin_data, output_bin_data

Nmlt   = 0
iFname = []
opera  = []
oFname = None
factor = 1

###### Functions for arguments
def usage(ARGV):
    print
    print "usage: %s [ofname] [factor] [x or /] [ifname1] [x or /] [ifname2] ..." % ARGV[0]
    print
    quit()

def set_args(ARGC, ARGV):
    global iFname, opera, oFname, Nmlt, factor
    
    oFname = ARGV[1].strip()
    factor = float(ARGV[2].strip())
    Nmlt   = (ARGC - 3) / 2
    
    for imlt in range(3, ARGC, 2):
        opera.append( ARGV[imlt+0].strip())
        iFname.append(ARGV[imlt+1].strip())
    
    for imlt in range(Nmlt):
        if (opera[imlt] != 'x' and opera[imlt] != '/' ):
            print("\nERROR: Invalid operator '%s', exit.\n" % opera[imlt]); quit()

###### Main part
if __name__ == "__main__":
    from numpy import empty
    
    argv = sys.argv; argc = len(argv)
    
    if (argc < 5 or (argc-1) % 2 != 0):
        usage(argv)
    
    set_args(argc, argv)#; print oFname; print factor; print Nmlt; print opera; print iFname; quit()
    
    yData, xData, Dummy = input_bin_data(iFname[0])
    if (xData is yData is Dummy is None):
        quit()
    
    Nconf = len(yData[0, :])
    Ndata = len(yData[:, 0])
    del yData; yData = empty((Ndata, Nconf))
    
    for idata in range(Ndata):
        for iconf in range(Nconf):
            yData[idata, iconf] = factor
    
    for imlt in range(Nmlt):
        tmpyData, tmpxData, Dummy = input_bin_data(iFname[imlt])
        if (tmpxData is tmpyData is Dummy is None):
            quit()
        
        for idata in range(Ndata):
            if (xData[idata] != tmpxData[idata]):
                print("\nERROR: xData is differ, exit.\n"); quit()
            for iconf in range(Nconf):
                if   (opera[imlt] == "x"):
                    yData[idata, iconf] *= tmpyData[idata, iconf]
                elif (opera[imlt] == "/"):
                    yData[idata, iconf] /= tmpyData[idata, iconf]
                else:
                    print("\nERROR: Invalid operator '%s', exit.\n" % opera[imlt]); quit()
        
        del tmpyData, tmpxData
    
    output_bin_data(oFname, yData, xData)
