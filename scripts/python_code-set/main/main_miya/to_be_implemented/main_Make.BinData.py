#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np
from MiscFuncs.DataStatistics import make_mean_err
from MiscFuncs.IO_TextData    import input_text_data
from MiscFuncs.IO_BinData     import output_bin_data

iFname = None
oFname = None

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [ofname (bin-data)] [ifname (2-column text-data)] ...\n" % ARGV[0])
    quit()

def set_args(ARGC, ARGV):
    global iFname, oFname
    
    oFname = ARGV[1].strip()
    
    iFname = [None] * (ARGC-2)
    for iarg in range(ARGC-2):
        iFname[iarg] = ARGV[iarg+2].strip()

###### Main part
if __name__ == "__main__":
    from struct import pack
    
    argv = sys.argv; argc = len(argv)
    
    if (argc < 3):
        usage(argv)
    
    set_args(argc, argv)
    
    tmpData = input_text_data(iFname[0])
    if (len(tmpData[0, :]) != 2):
        print("ERROR: ONLY 2-column text data can be read, exit.")
        quit()
    
    Nfile = len(iFname)
    print("#.ifile = %d" % Nfile)
    
    Ndata = len(tmpData[:, 0])
    print("#.data  = %d" % Ndata)
    
    xData = np.empty(Ndata)
    eData = np.empty(Ndata)
    yData = np.empty((Ndata, Nfile))
    
    for idata in range(Ndata):
        xData[idata]    = tmpData[idata, 0]
        yData[idata, 0] = tmpData[idata, 1]
    
    for iarg in range(1, Nfile):
        print("Reading: %s" % iFname[iarg])
        tmpData = input_text_data(iFname[iarg])
        
        if (len(tmpData[0, :]) != 2):
            print("ERROR: ONLY 2-column text data can be read, exit.")
            quit()
        
        for idata in range(Ndata):
            yData[idata, iarg] = tmpData[idata, 1]
            
            if (xData[idata] != tmpData[idata, 0]):
                print("ERROR: x-data are different, exit.")
                quit()
    
    for idata in range(Ndata):
        eData[idata] = make_mean_err(yData[idata, :])[1]
    
    print("Writing: %s" % oFname)
    output_bin_data(oFname, yData, xData, eData)
