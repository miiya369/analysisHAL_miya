#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Takaya Miyamoto
# Brief : Add fit parameter files for gaussian potential
# Date  : Thu May 17 16:19:40 JST 2018

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../..")
import numpy as np
import time

from fitting.io_params import input_params, output_params

### ================== Global Parameters Init. ================= ###
ofname  = None
ifnames = ["Pot.ScN12__ScN12_.V_C.J_0.t13.bin_size57.fitparam", 
           "Pot.ScN12__ScN12_.eff.3S1.t13.bin_size57.fitparam", 
           "Pot.ScN32__single.V_C.J_0.t13.bin_size57.fitparam", 
           "Pot.ScN32__single.eff.3S1.t13.bin_size57.fitparam"]
coeffs  = np.array([1.0/12.0, 3.0/12.0, 2.0/12.0, 6.0/12.0])
### =========================== Main =========================== ###
def main():
    Nadd = len(ifnames)
    if (Nadd == 0):
        print("\nERROR: No input files, exit."); return -1
    
    func_names =          [input_params(ifnames[i])[0] for i in range(Nadd)]
    iparams    = np.array([input_params(ifnames[i])[1] for i in range(Nadd)])
    Ngauss     = np.array([len(iparams[i,0,:]) // 2    for i in range(Nadd)])
    
    Nconf = len(iparams[0,:,0])
    
    for i in range(Nadd):
        for c in range(Nconf):
            for k in range(Ngauss[i]):
                iparams[i,c,2*k] *= coeffs[i]
    
    oparam = iparams[0,:,:]
    for i in range(1,Nadd):
        oparam  = np.append(oparam, iparams[i,:,:], axis=1)
    
    output_params(ofname, "12G", oparam)
    
    return 0

### ============================================================ ###
###### Functions for arguments ######
def usage(ARGV0):
    print("usage  : python %s [ofname] {options}" % os.path.basename(ARGV0))
    print("options:")
    print("      none")
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ofname =", ofname)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ofname
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ofname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--none'):
                pass
            elif (ARGV[i] == '--none'):
                pass
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
    
    if (main() != 0):
        exit("ERROR EXIT.")
