#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname = None
### =========================== Main =========================== ###

def main():
    from common.statistics         import make_mean_err
    from Tmatrix.io_Tmatrix        import input_Tmatrix
    from Tmatrix.convert_mat       import convert_TtoS
    from Tmatrix.calc_phase_shift  import calc_phase_Sii, within_one
    
### Input the T-matrix ###
    xData, Tmat = input_Tmatrix(ifname)
    if (xData is Tmat is None):
        return -1
    Nconf = len(Tmat[:,0,0,0])
    Ndata = len(Tmat[0,:,0,0])
    Nch   = len(Tmat[0,0,:,0])
    
### Convert T-matrix to S-matrix ###
    Smat = np.array([[convert_TtoS(Tmat[iconf,idata,:,:]) for idata in range(Ndata)] for iconf in range(Nconf)])
    
    if (False): # For Debug (Unitarity check for S-matrix)
        for iconf in range(Nconf):
            for idata in range(Ndata):
                Mat = np.dot(Smat[iconf,idata,:,:], np.conjugate(Smat[iconf,idata,:,:].T))
                print("%lf:\n" % xData[idata], Mat, "\n")
        return 0
    
### Phase shfit calculation for each Nch ###
    if   (Nch == 1):
        phase = np.array([calc_phase_Sii(Smat[i,:,0,0]) for i in range(Nconf)])
        print("# Plot index: Energy(MeV)|pha-m|phs-e|Abs(S)-m|Abs(S)-e")
        for idata in range(Ndata):
            print("%lf %e %e %e %e" % 
                  (xData[idata], *make_mean_err(phase[:,0,idata]), *make_mean_err(phase[:,1,idata])))
    elif (Nch == 2):
        phase  = np.array([[calc_phase_Sii(Smat[i,:,ich,ich], shift_disc = 50) for i in range(Nconf)]
                           for ich in range(Nch)])
        
        MixAng = np.array([[np.arccos(within_one(phase[0,i,1,iE])) * 90.0 / np.pi for iE in range(Ndata)] 
                           for i in range(Nconf)])
        
        phase[1,:,0,:] = np.array([[0.0 if (abs(phase[0,i,1,iE]-1.0) < 1e-10) else phase[1,i,0,iE] 
                                    for iE in range(Ndata)] 
                                   for i in range(Nconf)])
        print("# Plot index: Energy(MeV)|pha1-m|phs1-e|phs2-m|phs2-e|MixAng-m|MixAng-e|Inela-m|Inela-e")
        for idata in range(Ndata):
            print("%lf %e %e %e %e %e %e %e %e" % 
                  (xData[idata],
                   *make_mean_err(phase[0,:,0,idata]), *make_mean_err(phase[1,:,0,idata]),
                   *make_mean_err(MixAng[:,idata])   , *make_mean_err(phase[0,:,1,idata])))
    else:
        print("\nN.channel = %d has not been implemented yet, exit.\n" % Nch)
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    exit("usage: python %s [ifile]" % os.path.basename(ARGV0))

def check_args():
    print("# === Check Arguments ===")
    print("# ifile =", ifname)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname
    
    if (ARGC != 2):
        usage(ARGV[0])

    ifname = ARGV[1].strip()
    
    check_args()

### ============================================================ ###
### ============================================================ ###
if (__name__ == "__main__"):
    argv = sys.argv; argc = len(argv)
    if (argc == 1):
        usage(argv[0])
    
    set_args(argc, argv)
    
    t_start = time.time()
    if (main() != 0):
        exit("ERROR EXIT.")
    print("#\n# Elapsed time [s] = %d" % (time.time() - t_start))
