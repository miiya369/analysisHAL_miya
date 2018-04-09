#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np

from common.misc          import frange
from common.statistics    import make_mean_err
from fitting.io_params    import input_params
from fitting.fitfunc_type import set_fitfunc_from_fname

### ================== Global Parameters Init. ================= ###
ifname = None
r_min  = 0.001
r_del  = 0.01
r_max  = 2.5
### =========================== Main =========================== ###
def main():
    Fname, Params = input_params(ifname)
    if (Fname is Params is None):
        return -1
    
    Nconf   = len(Params[:,0])
    Nparam  = len(Params[0,:])
    FitFunc = set_fitfunc_from_fname(Fname)
    
    for r in frange(r_min, r_max, r_del):
        print("%lf %1.16e %1.16e" % 
              (r, *make_mean_err(np.array([FitFunc(r,*Params[iconf,:])
                                           for iconf in range(Nconf)]))))
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : %s [ifile] {options}" % os.path.basename(ARGV0))
    print("options:")
    print("      --r_min [Minimum range  (fm)] Default =", r_min)
    print("      --r_del [Range division (fm)] Default =", r_del)
    print("      --r_max [Maximum range  (fm)] Default =", r_max)
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile =", ifname)
    print("# r min =", r_min)
    print("# r del =", r_del)
    print("# r max =", r_max)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname, r_min, r_del, r_max
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--r_min'):
                r_min = float(ARGV[i+1])
            elif (ARGV[i] == '--r_del'):
                r_del = float(ARGV[i+1])
            elif (ARGV[i] == '--r_max'):
                r_max = float(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV[0])
    
    check_args()

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)

    if (argc == 1):
        usage(argv[0])

    set_args(argc, argv)

    if (main() != 0):
        exit("ERROR EXIT.")
