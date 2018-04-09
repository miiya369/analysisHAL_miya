#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Wed Feb  8 16:21:19 JST 2017
### Brief : Solve the schrodinger equation by using the Gauss Expansion method.

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from Fitting.IO_Params        import input_params
from SchGaussExp.PrintResults import print_results
from SchGaussExp.SolveSch     import solve_sch_GEM, calc_from_HamMat
from SchGaussExp.IO_HamMat    import input_HamMat

Fname   = None
mass    = 500.0
max_r   = 10.0
range_a = 0.8
Nbase   = 10
Nret    = 1
Np      = 0

### For Single-Folidng Potentials (2pF-type (WS-type) density)
from FoldingPotential.DensityFunc import dens_woods_saxon, calc_rho0

SFP_flg = True
#SFP_flg = False
tmpNA   = 208
tmpRA   = 6.80
tmpaA   = 0.515

tmpFactor = 1.0 # For test
###

Do_stochastic = False
EigValOnly    = True
FromHamMat    = False

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [Fit Param] {--option}\n" % ARGV[0])
    print("option:")
    print("      --mass        [reduced mass          (MeV)] Default ="),; print mass
    print("      --Np          [#.charge                   ] Default ="),; print Np
    print("      --max_r       [Max range of gaussian ( fm)] Default ="),; print max_r
    print("      --range_a     [The base of range          ] Default ="),; print range_a
    print("      --Nbase       [#.base                     ] Default ="),; print Nbase
    print("      --Nret        [#.eigen value to calculate ] Default ="),; print Nret
    print("      --stochastic   Using stochastic variational method"),; print
    print("      --out_wave     Output wave function"),; print
    print("      --from_mat     Input Hamiltonian matrixes"),; print
    print; quit()

def set_args(ARGC, ARGV):
    global Fname, mass, max_r, range_a, Nbase, Nret, Do_stochastic, EigValOnly, FromHamMat
    global tmpNA, tmpRA, tmpaA, Np, tmpFactor
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    Fname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass'):
                mass = float(ARGV[i+1])
            elif (ARGV[i] == '--max_r'):
                max_r = float(ARGV[i+1])
            elif (ARGV[i] == '--range_a'):
                range_a = float(ARGV[i+1])
            elif (ARGV[i] == '--Nbase'):
                Nbase = int(ARGV[i+1])
            elif (ARGV[i] == '--Nret'):
                Nret = int(ARGV[i+1])
            elif (ARGV[i] == '--A'):
                tmpNA = int(ARGV[i+1])
            elif (ARGV[i] == '--RA'):
                tmpRA = float(ARGV[i+1])
            elif (ARGV[i] == '--aA'):
                tmpaA = float(ARGV[i+1])
            elif (ARGV[i] == '--Np'):
                Np = int(ARGV[i+1])
            elif (ARGV[i] == '--out_wave'):
                EigValOnly = False
            elif (ARGV[i] == '--stochastic'):
                Do_stochastic = True
            elif (ARGV[i] == '--from_mat'):
                FromHamMat = True
            elif (ARGV[i] == '--factor'):
                tmpFactor = float(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     ="),; print Fname
    print("# mass      ="),; print mass
    print("# Np        ="),; print Np
    print("# max range ="),; print max_r
    print("# range a   ="),; print range_a
    print("# N.base    ="),; print Nbase
    print("# N.retarn  ="),; print Nret
    print("# Do_stoc.  ="),; print Do_stochastic
    print("# Out wave  ="),; print (not EigValOnly)
    print("# From mat  ="),; print FromHamMat
    print("# =======================")

###### Main part
if __name__ == "__main__":
    import numpy as np
    import re
    import time; JOB_start = time.time()
    
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv)
    
    if (FromHamMat):
        mass          = "From mat"
        max_r         = "From mat"
        range_a       = "From mat"
        Nbase         = "From mat"
        Do_stochastic = "From mat"
    
    check_args()
    
    if (SFP_flg):
        tmpRho0 = calc_rho0(tmpNA, dens_woods_saxon, np.array((tmpRA, tmpaA)))
        print("# NA        ="),; print tmpNA
        print("# RA        ="),; print tmpRA
        print("# aA        ="),; print tmpaA
        print("# rho0      ="),; print tmpRho0
        if (tmpFactor != 1): # For Test
            print("# factor    ="),; print tmpFactor
    
### Input parameters ###
    if (Fname != 'DEBUG'):
        if (FromHamMat):
            range_a, max_r, Range, Ham = input_HamMat(Fname)
            if (Range is Ham is None):
                quit()
        
        else:
            Ftype, Params = input_params(Fname)
            if (Ftype is Params is None):
                quit()
                
            Nconf = len(Params[0, :])
    
    else:
        ### For Debug
        Vzero = -50.0; Rzero = 2.0; mass = 500.0; Debug_Expect = (-9.85391)
        #Vzero = -100.0; Rzero = 4.0; mass = 1000.0; Debug_Expect = (-90.291, -61.6334, -16.6833)
        print("# DEGUB MODE...")
        print("# V_0 (MeV) ="),; print Vzero
        print("# R_0 ( fm) ="),; print Rzero
        print("# Func type = 3-dim Idogata")
        print("# Expect En ="),; print Debug_Expect
        SFP_flg   = False
        Ftype     = "SW"
        Nconf     = 1
        Params    = np.empty((2, Nconf))
        Params[0] = Vzero; Params[1] = Rzero
    
    if (FromHamMat):
        if (not SFP_flg):
            Eigen_val, Eigen_vec = calc_from_HamMat(Nret, Range, Ham, Np, EigValOnly)
        else:
            Eigen_val, Eigen_vec = calc_from_HamMat(Nret, Range, Ham,
                                                    np.array((Np, tmpRho0, tmpRA, tmpaA, tmpNA)),
                                                    EigValOnly)
        if (Eigen_val is Eigen_vec is None):
            quit()
        
        print_results(Eigen_val, Eigen_vec, Range, max_r, 0.01, EigValOnly)
        
        print("#\n# JOB END: elapsed_time = %lf [sec]" % (time.time() - JOB_start))
        quit()
    
    if (Do_stochastic):
### The stocastic variational method ###
        Range = np.empty(0); tmp_val = 1e+30
        for i in range(Nbase):
            Range = np.append(Range, 0.0)
            while(True):
                Range[i] = np.random.random() * max_r; #print Range
                
                if (not SFP_flg):
                    val = solve_sch_GEM(1, Range, Params[:, 0], mass, Ftype, Np, True)[0]
                else:
                    val = solve_sch_GEM(1, Range, Params[:, 0],
                                        np.array((mass, tmpRho0, tmpRA, tmpaA, tmpNA)),
                                        Ftype, Np, True)[0]
                if (val is None):
                    quit()
                
                if (val[0] < tmp_val):
                    tmp_val = val[0]
                    break
    else:
### Using the Nbase-Gaussian (Range is set as geometric series) ###
        Range = np.empty(Nbase)
        for i in range(Nbase):
            Range[i] = max_r * pow(range_a, i)
    
    Eigen_val = np.empty((       Nret, Nconf))
    Eigen_vec = np.empty((Nbase, Nret, Nconf))
    
    print("#")
    for i in range(Nconf):
        if (not SFP_flg):
            val, vec = solve_sch_GEM(Nret, Range, Params[:, i], mass, Ftype, Np, EigValOnly)
        else:
            val, vec = solve_sch_GEM(Nret, Range, Params[:, i], 
                                     np.array((mass, tmpRho0, tmpRA, tmpaA, tmpNA)),
                                     Ftype, Np, EigValOnly)
        
        if   (val is vec is None):
            quit()
        elif (EigValOnly):
            for iret in range(Nret):
                Eigen_val[iret, i] = val[iret]
        else:
            for iret in range(Nret):
                Eigen_val[iret, i] = val[iret]
                for ibase in range(Nbase):
                    Eigen_vec[ibase, iret, i] = vec[ibase, iret]
        
        print("# conf = %d: End" % i)
    
    print_results(Eigen_val, Eigen_vec, Range, max_r, 0.01, EigValOnly)
    
    print("#\n# JOB END: elapsed_time = %lf [sec]" % (time.time() - JOB_start))
