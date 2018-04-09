#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Wed Feb  8 16:21:19 JST 2017
### Brief : Solve Schrodinger equation with 1st perturbation
###       : For the folding potential of Coulomb + Fourier-Bessel Coefficient

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from SchPert.IO_wave          import input_wave
from SchPert.SolveSch_Pert    import solve_sch_1st_pert_Coulomb_FBcoe

iFname_wave = None

NA = 208

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [input file (Output of Sch.Gauss.Exp)] {--option}\n" % ARGV[0])
    print("option:")
    print("      --A [#.nucleus] Default ="),; print NA
    print; quit()

def set_args(ARGC, ARGV):
    global iFname_wave, NA
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFname_wave = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--A'):
                NA = int(ARGV[i+1])
            elif (ARGV[i] == '--A'):
                NA = int(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile wave ="),; print iFname_wave
    print("# N.nucleus  ="),; print NA
    print("# =======================")

###### Main part
if __name__ == "__main__":
    import numpy as np
    import time; JOB_start = time.time()
        
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv); check_args()
    
### Input Wave function params
    Range, Coeff = input_wave(iFname_wave);
    if (Range is Coeff is None):
        quit()
    
    Nconf = len(Coeff[0, :]) #; print Range; print Coeff; quit()
        
### Calculation & Output    
    Energy_1 = np.empty(Nconf)
    print("#")
    for i in range(Nconf):
        Energy_1[i] = solve_sch_1st_pert_Coulomb_FBcoe(Range, Coeff[:, i], NA)
        print("# conf = %d: End" % i)
    
    Mean, Err = make_mean_err(Energy_1)
    print("#\n# E_1 = %1.16f +/- %1.16f" % (Mean, Err))
    
    print("#\n# JOB END: elapsed_time = %lf [sec]" % (time.time() - JOB_start))
