#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np

from common.io_data_bin import input_bin_data, output_bin_data

### ================== Global Parameters Init. ================= ###

### =========================== Main =========================== ###
def main(ofpath, ifpaths):
    Naverage = len(ifpaths)
    
    print("#.average =", Naverage); print("")
    
    y_data0, x_data0, dummy = input_bin_data(ifpaths[0])
    if (y_data0 is x_data0 is dummy is None):
        return -1
    
    for i in range(1, Naverage):
        y_datai, x_datai, dummy = input_bin_data(ifpaths[i])
        if (y_datai is x_datai is dummy is None):
            return -1
        
        if ((x_data0 != x_datai).any()):
            print("ERROR: x_data in '%s' is differ, exit." % ifpaths[i])
            return -1
        
        y_data0 += y_datai
    
    y_data0 /= float(Naverage)
    
    print(""); output_bin_data(ofpath, y_data0, x_data0)
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    exit("usage  : %s [ofname] [ifname1] [ifname2] ..." % os.path.basename(ARGV0))

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc < 3):
        usage(argv[0])
    
    if (main(argv[1], argv[2:]) != 0):
        exit("ERROR EXIT.")
