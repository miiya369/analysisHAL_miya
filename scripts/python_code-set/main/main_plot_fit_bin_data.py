#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np

from common.io_data_bin import input_bin_data

### ================== Global Parameters Init. ================= ###

### =========================== Main =========================== ###
def main(ifpath):
    y_data, x_data, e_data = input_bin_data(ifpath)
    if (y_data is x_data is e_data is None):
        return -1
    
    for i in range(len(x_data)):
        y_data_m = np.average(y_data[:,i])
        print(x_data[i], y_data_m, e_data[i])
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    exit("usage  : %s [ifname]" % os.path.basename(ARGV0))

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc != 2):
        usage(argv[0])
    
    if (main(argv[1]) != 0):
        exit("ERROR EXIT.")
