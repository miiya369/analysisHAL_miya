#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np
import time

from common.io_data_bin            import input_bin_data_nxyz
from sch_finiteV.solve_sch_finiteV import solve_sch_Fvol
from sch_finiteV.print_results     import print_results, output_results
from lattice.phase_shift_Luscher   import k_cot_d_Luscher
from misc_QM.scattering_Idogata    import tan_d_Idogata3D, Eb_Idogata3D_Swave

### ================== Global Parameters Init. ================= ###
ifname  = None
rmass   = 0.6215
lat_a   = None
Nret    = 40
Lsize   = 16
ofbase  = None

EigValOnly = True

hbar_c  = 197.327053

### For Debug
Vzero   = -200.0
Rzero   = 1.0
mu      = 500
dummy_d = 0.0
### =========================== Main =========================== ###
def main():
    t_start = time.time()
    
### Input parameters ###
    if (ifname != 'DEBUG'):
        Pot = input_bin_data_nxyz(ifname, data_type = "float")
        if (Pot is None):
            return -1
        Nconf = len(Pot[:,0,0,0])
    else:
        ### For Debug
        print("# DEGUB MODE...")
        print("# V_0 (MeV) =", Vzero)
        print("# R_0 ( fm) =", Rzero)
        print("#  mu (MeV) =",    mu)
        print("#   a ( fm) =", lat_a)
        print("# Func type = 3-dim Idogata")
        print("# Expect En =", Eb_Idogata3D_Swave(Vzero, Rzero, mu))
        #test(Vzero, Rzero, mu, dummy_d * lat_a / hbar_c); return 0
        Nconf  = 1
        Pot    = np.zeros((Nconf, Lsize, Lsize, Lsize))
        for i in range(Nconf):
            for z in range(Lsize):
                for y in range(Lsize):
                    for x in range(Lsize):       
                        if (np.sqrt(x**2 + y**2 + z**2) < Rzero/lat_a):
                            Pot[i,z,y,x] = Vzero * lat_a / hbar_c
    
    Results = np.array([])
    print("#")
    for i in range(Nconf):
        Results = np.append(Results, solve_sch_Fvol(Pot[i,:,:,:],
                                                    rmass, Lsize, Nret, not EigValOnly))
        print("# conf = %03d/%03d end." % (i+1, Nconf))
    
    if (EigValOnly):
        EigVal = np.reshape(Results, (Nconf, Nret))
        EigVec = None
    else:
        EigTmp = np.reshape(Results, (Nconf,    2))
        EigVal = np.array([EigTmp[i,0] for i in range(Nconf)])
        EigVec = np.array([EigTmp[i,1] for i in range(Nconf)])
        
        for i in range(Nconf):
            for k in range(Nret):
                if (EigVec[i,1,k] < 0):
                    EigVec[i,:,k] *= -1.0
    
    print_results(EigVal, EigVec, lat_a, EigValOnly)
    print("#")
    
    if (ofbase is not None):
        output_results(ofbase, EigVal, EigVec, lat_a, EigValOnly)
    
    print("#\n# Elapsed time [s] = %d" % (time.time() - t_start))
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --rmass       [reduced mass (Lattice Unit)] Default =", rmass)
    print("      --lat_a       [Lattice spacing        (fm)] Default =", lat_a)
    print("      --Nret        [#.eigen value to calculate ] Default =", Nret)
    print("      --Lsize       [Lattice size  to calculate ] Default =", Lsize)
    print("      --out_wave     Output wave function")
    print("      --ofbase      [Base name for output files ] Default =", ofbase)
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     =", ifname)
    print("# rmass     =", rmass)
    print("# lat.space =", lat_a)
    print("# Lsize     =", Lsize)
    print("# N.retarn  =", Nret)
    print("# out wave  =", (not EigValOnly))
    print("# ofbase    =", ofbase)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname, rmass, lat_a, Nret, Lsize, EigValOnly, dummy_d, ofbase
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--rmass'):
                rmass = float(ARGV[i+1])
            elif (ARGV[i] == '--lat_a'):
                lat_a = float(ARGV[i+1])
            elif (ARGV[i] == '--Lsize'):
                Lsize = int(ARGV[i+1])
            elif (ARGV[i] == '--Nret'):
                Nret = int(ARGV[i+1])
            elif (ARGV[i] == '--out_wave'):
                EigValOnly = False
            elif (ARGV[i] == '--ofbase'):
                ofbase = ARGV[i+1].strip()
            elif (ARGV[i] == '--test'):
                dummy_d = float(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV[0])
    
    if (ifname == 'DEBUG'):
        rmass = mu * lat_a / hbar_c
    
    check_args()

def test(V0, R0, Mu, EigVal):
    #if (True):
    if (False):
        for k in range(0, 1000, 1):
            print(k**2,
                  k/tan_d_Idogata3D(k, V0, R0, Mu, 0),
                  k_cot_d_Luscher((k * lat_a / hbar_c)**2, Lsize) * hbar_c / lat_a)
    else:
        Mu *= lat_a / hbar_c
        k2  = 2.0*Mu*EigVal * (hbar_c / lat_a)**2
        print(k2, k_cot_d_Luscher(2.0*Mu*EigVal, Lsize) * hbar_c / lat_a)

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv[0])
    
    set_args(argc, argv)
    
    if (main() != 0):
        exit("ERROR EXIT.")
