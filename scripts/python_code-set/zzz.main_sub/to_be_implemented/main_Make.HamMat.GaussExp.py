#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Wed Feb  8 16:21:19 JST 2017
### Brief : Make Hamiltonian matrixes for the Gauss Expansion method.

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from Fitting.IO_Params        import input_params
from SchGaussExp.SetMatrix    import *
from SchGaussExp.IO_HamMat    import input_HamMat, output_HamMat

iFname  = None
oFname  = "/dev/null"
mass    = 500.0
max_r   = 10.0
range_a = 0.8
Nbase   = 10

### For Single-Folidng Potentials
from FoldingPotential.DensityFunc import dens_woods_saxon, calc_rho0

SFP_flg = True
#SFP_flg = False
tmpNA   = 208
tmpRA   = 6.80
tmpaA   = 0.515
###

Do_stochastic = False
FromHamMat    = False

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [Fit params] {--option}\n" % ARGV[0])
    print("option:")
    print("      --ofile       [output file name           ] Default ="),; print oFname
    print("      --mass        [reduced mass          (MeV)] Default ="),; print mass
    print("      --max_r       [Max range of gaussian ( fm)] Default ="),; print max_r
    print("      --range_a     [The base of range          ] Default ="),; print range_a
    print("      --Nbase       [#.base                     ] Default ="),; print Nbase
    print("      --stochastic   Using stochastic variational method"),; print
    print("      --from_mat     Input Hamiltonian matrixes"),; print
    print; quit()

def set_args(ARGC, ARGV):
    global iFname, oFname, mass, max_r, range_a, Nbase, Do_stochastic, FromHamMat
    global tmpNA, tmpRA, tmpaA
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass'):
                mass = float(ARGV[i+1])
            elif (ARGV[i] == '--ofile'):
                oFname = ARGV[i+1].strip()
            elif (ARGV[i] == '--max_r'):
                max_r = float(ARGV[i+1])
            elif (ARGV[i] == '--range_a'):
                range_a = float(ARGV[i+1])
            elif (ARGV[i] == '--Nbase'):
                Nbase = int(ARGV[i+1])
            elif (ARGV[i] == '--A'):
                tmpNA = int(ARGV[i+1])
            elif (ARGV[i] == '--RA'):
                tmpRA = float(ARGV[i+1])
            elif (ARGV[i] == '--aA'):
                tmpaA = float(ARGV[i+1])
            elif (ARGV[i] == '--stochastic'):
                Do_stochastic = True
            elif (ARGV[i] == '--from_mat'):
                FromHamMat = True
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     ="),; print iFname
    print("# ofile     ="),; print oFname
    print("# mass      ="),; print mass
    print("# max range ="),; print max_r
    print("# range a   ="),; print range_a
    print("# N.base    ="),; print Nbase
    print("# Do_stoc.  ="),; print Do_stochastic
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
    check_args()
    
    if (SFP_flg):
        tmpRho0 = calc_rho0(tmpNA, dens_woods_saxon, np.array((tmpRA, tmpaA)))
        print("# NA        ="),; print tmpNA
        print("# RA        ="),; print tmpRA
        print("# aA        ="),; print tmpaA
        print("# rho0      ="),; print tmpRho0
    
### Input parameters ###
    if (FromHamMat):
        print(" ERROR: Have not implemented yet, exit."); quit()
    else:
        Ftype, Params = input_params(iFname)
        if (Ftype is Params is None):
            quit()
        
        Nconf = len(Params[0, :])
    
    if (Do_stochastic):
### The stocastic variational method ###
        Range = np.empty(0); tmp_val = 1e+30
        for i in range(Nbase):
            Range = np.append(Range, 0.0)
            while(True):
                Range[i] = np.random.random() * max_r; #print Range
                
                if (not SFP_flg):
                    val = solve_sch_GEM(1, Range, Params[:, 0], mass, Ftype, 0, True)[0]
                else:
                    val = solve_sch_GEM(1, Range, Params[:, 0],
                                        np.array((mass, tmpRho0, tmpRA, tmpaA, tmpNA)),
                                        Ftype, 0, True)[0]
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
    
    HamMat = np.empty((Nbase, Nbase, Nconf))
    
    print("#")
    for iconf in range(Nconf):
        if (not SFP_flg):
            if   (re.match('^[1-9]G$', Ftype) is not None):
                tmp_Ham = set_mat_gauss(Range, Params[:, iconf], mass)[0]
                #tmp_Ham = set_mat(Range, Params[:, iconf], Ftype, mass)[0]
            
            else:
                tmp_Ham = set_mat(Range, Params[:, iconf], Ftype, mass)[0]
                
        else:
            if   (re.match('^[1-9]G$', Ftype) is not None):
                tmp_Ham = set_mat_gauss_fpot_WS(Range, Params[:, iconf], mass, 
                                                np.array((tmpRho0, tmpRA, tmpaA)))[0]
                #tmp_Ham = set_mat_fpot_WS(Range, Params[:, iconf], Ftype, mass, 
                #                          np.array((tmpRho0, tmpRA, tmpaA)))[0]
            
            else:
                tmp_Ham = set_mat_fpot_WS(Range, Params[:, iconf], Ftype, mass, 
                                          np.array((tmpRho0, tmpRA, tmpaA)))[0]
        
        for i in range(Nbase):
            for j in range(Nbase):
                HamMat[i, j, iconf] = tmp_Ham[i, j]
    
        print("# conf = %d: End" % iconf)
    
    output_HamMat(oFname, range_a, max_r, Range, HamMat)
    
    print("#\n# JOB END: elapsed_time = %lf [sec]" % (time.time() - JOB_start))
