#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname  = None
mass    = 0.5
lat_a   = None
Nstat   = 3
Lsize   = 16
Nproc   = 1

eval_only = True
ofbase    = None

hbar_c  = 197.327053
dummy_d = None

def fproc_wrapper(args):
    from sch_finiteV.solve_sch_finiteV import solve_sch_Fvol
    a_iconf     = args[0]
    a_pot       = args[1]
    a_mass      = args[2]
    a_Lsize     = args[3]
    a_Nstat     = args[4]
    a_eval_only = args[5]
    results = solve_sch_Fvol(a_pot, a_mass, a_Lsize, a_Nstat, not a_eval_only)
    print("# Solve Schrodinger equation in the finite volume... %3d end" % a_iconf)
    return results
### =========================== Main =========================== ###
def main():
    from common.io_data_bin      import input_bin_data_nxyz
    from sch_finiteV.print_eigen import print_eigen, output_eigen
    
### Input parameters ###
    if (ifname != 'DEBUG'):
        pot = input_bin_data_nxyz(ifname, data_type = "float")
        if (pot is None):
            return -1
        Nconf = len(pot[:,0,0,0])
    else:
        ### For Debug
        from misc_QM.scattering_Idogata import Eb_Idogata3D_Swave
        global mass, lat_a
        if (lat_a is None):
            lat_a = 0.5
        Vzero = -50.0; Rzero = 2.0; mass = 500.0
        #Vzero = -100.0; Rzero = 4.0; mass = 1000.0
        Nconf = 1
        pot   = np.array([[[[Vzero * lat_a/hbar_c if (np.sqrt(x**2+y**2+z**2) < Rzero/lat_a) else 0.0
                             for x in range(Lsize)] for y in range(Lsize)] for z in range(Lsize)]
                          for i in range(Nconf)])
        print("# DEGUB MODE...")
        print("# V_0 (MeV) =", Vzero)
        print("# R_0 ( fm) =", Rzero)
        print("# mass(MeV) =", mass )
        print("#   a ( fm) =", lat_a)
        print("# Func type = 3-dim Idogata")
        print("# Expect En =", Eb_Idogata3D_Swave(Vzero, Rzero, mass))
        if (dummy_d is not None):
            test(Vzero, Rzero, mass, dummy_d * lat_a/hbar_c); return 0
        else:
            mass *= lat_a/hbar_c
    
### Solve the Schrodinger equation in the finite volume ###
    print("#\n# Solve Schrodinger equation in the finite volume...")
    if (Nproc == 1):
        eigens = np.array([fproc_wrapper((iconf, pot[iconf], mass, Lsize, Nstat, eval_only)) for iconf in range(Nconf)])
    else:
        args_procs = [(iconf, pot[iconf], mass, Lsize, Nstat, eval_only) for iconf in range(Nconf)]
        with Pool(Nproc) as proc:
            eigens = np.array(proc.map(fproc_wrapper, args_procs))
    print("# Solve Schrodinger equation in the finite volume... all end")
    
### print results ###
    if (eval_only):
        eval = eigens
        evec = None
    else:
        eval = np.array([eigens[iconf,0] for iconf in range(Nconf)])
        evec = np.array([eigens[iconf,1] for iconf in range(Nconf)])
        for i in range(Nconf):
            for k in range(Nstat):
                if (evec[i,1,k] < 0):
                    evec[i,:,k] *= -1.0
    
    print_eigen(eval, evec, lat_a, eval_only)
    if (ofbase is not None):
        output_eigen(ofbase, eval, evec, lat_a, eval_only)
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --mass     [reduced mass (Lattice Unit)] Default =", mass)
    print("      --lat_a    [Lattice spacing        (fm)] Default =", lat_a)
    print("      --Nstat    [#.eigen state to calculate ] Default =", Nstat)
    print("      --Lsize    [Lattice size  to calculate ] Default =", Lsize)
    print("      --Nproc    [#.process                  ] Default =", Nproc)
    print("      --out_wave  Output wave function")
    print("      --ofbase   [Base name for output files ] Default =", ofbase)
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     =", ifname)
    print("# mass      =", mass)
    print("# lat.space =", lat_a)
    print("# Lsize     =", Lsize)
    print("# N.state   =", Nstat)
    print("# N.process =", Nproc)
    print("# out wave  =", (not eval_only))
    print("# ofbase    =", ofbase)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname, mass, lat_a, Nstat, Lsize, Nproc, eval_only, dummy_d, ofbase
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass'):
                mass = float(ARGV[i+1])
            elif (ARGV[i] == '--lat_a'):
                lat_a = float(ARGV[i+1])
            elif (ARGV[i] == '--Lsize'):
                Lsize = int(ARGV[i+1])
            elif (ARGV[i] == '--Nstat'):
                Nstat = int(ARGV[i+1])
            elif (ARGV[i] == '--Nproc'):
                Nproc = int(ARGV[i+1])
            elif (ARGV[i] == '--out_wave'):
                eval_only = False
            elif (ARGV[i] == '--ofbase'):
                ofbase = ARGV[i+1].strip()
            elif (ARGV[i] == '--test'):
                dummy_d = float(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
    check_args()

def test(V0, R0, Mu, eval):
    from misc_QM.scattering_Idogata  import tan_d_Idogata3D
    from lattice.phase_shift_Luscher import k_cot_d_Luscher
    print("#")
    if (False):
        for k in range(0, 1000, 1):
            print(k**2,
                  k/tan_d_Idogata3D(k, V0, R0, Mu, 0),
                  k_cot_d_Luscher((k * lat_a/hbar_c)**2, Lsize) * hbar_c/lat_a)
    else:
        Mu *= lat_a/hbar_c
        k2  = 2.0*Mu*eval * (hbar_c/lat_a)**2
        print(k2, k_cot_d_Luscher(2.0*Mu*eval, Lsize) * hbar_c/lat_a)

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
