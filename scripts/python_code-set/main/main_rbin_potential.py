#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname = None
ofname = "/dev/null"
#rbsize = 0.01
rbsize = 0.0845808
### =========================== Main =========================== ###

def main():
    from common.io_data_bin import input_bin_data, output_bin_data
    
### Input data ###
    iy, ix, ie = input_bin_data(ifname)
    if (iy is ix is ie is None):
        return -1
    
    # For avoiding NaN in the center of potential (e.g. The case of tensor force)
    #Nd = len(ix)
    #iy = iy[:,1:Nd-1]
    #ix = ix[  1:Nd-1]
    
    Nc = len(iy[:,0])
    Nd = len(iy[0,:])

    oy = np.empty((Nc, 0))
    ox = np.empty(0)

    oyerr = np.empty(0)
    oxerr = np.empty(0)
    
### Binding for r-coordinate ###
    for r in np.arange(0.0, np.max(ix), rbsize):
        by  = np.zeros(Nc)
        by2 = np.zeros(Nc)
        bx  = 0
        bx2 = 0
        nb  = 0
        for ir in range(Nd):
            if (r <= ix[ir] < r+rbsize):
                by  += iy[:,ir]
                by2 += iy[:,ir]**2
                bx  += ix[  ir]
                bx2 += ix[  ir]**2
                nb  += 1
        
        if (nb != 0):
            oy = np.append(oy, (by/float(nb)).reshape((Nc,1)), axis=1)
            ox = np.append(ox, (bx/float(nb)))
            
            if (nb == 1):
                oyerr = np.append(oyerr, 0.0)
                oxerr = np.append(oxerr, 0.0)
            else:
                oyerr = np.append(oyerr, np.max(np.sqrt(by2/float(nb)-(by/float(nb))**2)) / np.sqrt(nb-1))
                oxerr = np.append(oxerr,        np.sqrt(bx2/float(nb)-(bx/float(nb))**2)  / np.sqrt(nb-1))
    
### Output results ###
    for ir in range(len(ox)):
        print("%lf %e %e %e %e" % (ox[ir], np.mean(oy[:,ir]), np.std(oy[:,ir])*np.sqrt(Nc-1),
                                   oyerr[ir], oxerr[ir]))
    
    output_bin_data(ofname, oy, ox)
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --ofname [output file name        ] Default =", ofname)
    print("      --rbsize [binsize for r-coordinate] Default =", rbsize)
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile =", ifname)
    print("# ofile =", ofname)
    print("# bsize =", rbsize)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname, ofname, rbsize
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--ofname'):
                ofname = ARGV[i+1].strip()
            elif (ARGV[i] == '--rbsize'):
                rbsize = float(ARGV[i+1].strip())
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
