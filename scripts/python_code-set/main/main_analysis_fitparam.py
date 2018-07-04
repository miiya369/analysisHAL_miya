#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname = None
r_min  = 0.001
r_del  = 0.01
r_max  = 2.5
### =========================== Main =========================== ###

def main():
    from common.misc          import frange
    from common.statistics    import make_mean_err
    from fitting.io_params    import input_params
    from fitting.fitfunc_type import set_fitfunc_from_fname
    
    func_name, params = input_params(ifname)
    if (func_name is params is None):
        return -1
    
    Nconf    = len(params[:,0])
    Nparam   = len(params[0,:])
    fit_func = set_fitfunc_from_fname(func_name)
    
    for r in frange(r_min, r_max, r_del):
        print("%lf %1.16e %1.16e" %
              (r, *make_mean_err(np.array([fit_func(r,*params[iconf,:]) for iconf in range(Nconf)]))))
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile] {options}\n" % os.path.basename(ARGV0))
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
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
    check_args()
### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    if (argc == 1):
        usage(argv[0])
    
    set_args(argc, argv)
    
    t_start = time.time()
    if (main() != 0):
        exit("ERROR EXIT.")
    print("#\n# Elapsed time [s] = %d" % (time.time() - t_start))
