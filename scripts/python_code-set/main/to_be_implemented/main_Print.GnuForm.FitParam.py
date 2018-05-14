#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np

from MiscFuncs.Misc       import *
from Fitting.IO_Params    import input_params
from Fitting.PrintGnuForm import print_gnu

iFname = None

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [input bin-data] {--option}\n" % ARGV[0])
    print("option:")
    print("      None")
    print; quit()

def set_args(ARGC, ARGV):
    global iFname
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--pass'):
                pass
            elif (ARGV[i] == '--pass'):
                pass
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile ="),; print iFname
    print("# =======================")

###### Main part
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv)
    
    check_args()
    
### Input data ###
    Ftype, Params = input_params(iFname)
    if (Ftype is Params is None):
        quit()

    Nconf = len(Params[0, :])
        
### Print Gnu-Form ###
    for iconf in range(Nconf):
        print_gnu(Ftype, Params[:, iconf])
