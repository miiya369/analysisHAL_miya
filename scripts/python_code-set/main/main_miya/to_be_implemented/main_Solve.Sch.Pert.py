#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Wed Feb  8 16:21:19 JST 2017
### Brief : Solve Schrodinger equation with the perturbation

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from SchPert.IO_wave          import input_wave
from SchPert.SolveSch_Pert    import solve_sch_0th_pert, solve_sch_1st_pert
from Fitting.IO_Params        import input_params
from Fitting.FitFunctionType  import set_fitfunc_from_fname

iFname_wave = None
iFname_orgV = None
iFname_perV = 'Coulomb'

mass = 2000.0
Np   = 82

### For Single-Folidng Potentials (2pF-type (WS-type) density)
from FoldingPotential.DensityFunc   import dens_woods_saxon, calc_rho0
from FoldingPotential.Integrand     import integrand_1D_PotOpt_DensWS
from FoldingPotential.SolveIntegral import solve_Fint_1D
from Fitting.FitFunctionForm        import fitfunc_Coulomb, fitfunc_Coulomb_wF

SFP_flg = True
#SFP_flg = False
tmpNA   = 208
tmpRA   = 6.8
tmpaA   = 0.515
###

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [input file (Output of Sch.Gauss.Exp)] {--option}\n" % ARGV[0])
    print("option:")
    print("      --mass [          reduced mass (MeV)] Default ="),; print mass
    print("      --Np   [#.charge (For Coulomb force)] Default ="),; print Np
    print("      --orgV [param file for 0th-potential] Default ="),; print iFname_orgV
    print("      --perV [param file for 1st-potential] Default ="),; print iFname_perV
    print; quit()

def set_args(ARGC, ARGV):
    global iFname_wave, iFname_orgV, iFname_perV, mass, Np
    global tmpNA, tmpRA, tmpaA
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFname_wave = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass'):
                mass = float(ARGV[i+1])
            elif (ARGV[i] == '--orgV'):
                iFname_orgV = ARGV[i+1].strip()
            elif (ARGV[i] == '--perV'):
                iFname_perV = ARGV[i+1].strip()
            elif (ARGV[i] == '--A'):
                tmpNA = int(ARGV[i+1])
            elif (ARGV[i] == '--Np'):
                Np = int(ARGV[i+1])
            elif (ARGV[i] == '--RA'):
                tmpRA = float(ARGV[i+1])
            elif (ARGV[i] == '--aA'):
                tmpaA = float(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile wave ="),; print iFname_wave
    print("# ifile orgV ="),; print iFname_orgV
    print("# ifile perV ="),; print iFname_perV
    print("# mass       ="),; print mass
    print("# N.charge   ="),; print Np
    print("# =======================")

###### Main part
if __name__ == "__main__":
    import numpy as np
    import time; JOB_start = time.time()
        
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv); check_args()
    
    if (SFP_flg):
        Params_dens = np.array((tmpRA, tmpaA))
        tmpRho0     = calc_rho0(tmpNA, dens_woods_saxon, Params_dens)
        print("# NA        ="),; print tmpNA
        print("# RA        ="),; print tmpRA
        print("# aA        ="),; print tmpaA
        print("# rho0      ="),; print tmpRho0
    
### Input Wave function params
    Range, Coeff = input_wave(iFname_wave);
    if (Range is Coeff is None):
        quit()
    
    Nconf = len(Coeff[0, :]) #; print Range; print Coeff; quit()
    
### Input potential params (For original V(r))
    if (iFname_orgV is not None):
        factor_orgV = 1.0
        
        FuncName_orgV, Params_orgV = input_params(iFname_orgV)
        if (FuncName_orgV is Params_orgV is None):
            quit()
        if (Nconf != len(Params_orgV[0, :])):
            print("\nERROR: #.conf is differ between Wave[i] and org.Pot[i], exit.\n"); quit()
        
        Func_orgV = set_fitfunc_from_fname(FuncName_orgV)
        if (Func_orgV is None):
            quit()
        
        if (SFP_flg):
            factor_orgV = tmpRho0
            
            tmpParams_orgV = Params_orgV
            tmpFunc_orgV   = Func_orgV
            
            Func_orgV   = lambda r,Ag1,Ag2,Ag3,Ag4: factor_orgV * solve_Fint_1D(r,Ag1,Ag2,Ag3,Ag4)
            Params_orgV = np.array([[None]*Nconf]*4)
            
            for i in range(Nconf):
                Params_orgV[0, i] = integrand_1D_PotOpt_DensWS
                Params_orgV[1, i] = tmpFunc_orgV
                Params_orgV[2, i] = tmpParams_orgV[:, i]
                Params_orgV[3, i] = Params_dens
    
    #print factor_orgV; print Func_orgV; print Params_orgV; quit()
    
### Input potential params (For perturbative V(r))
    factor_perV = 1.0
    
    if (iFname_perV == 'Coulomb'):
        #Func_perV   = fitfunc_Coulomb
        Func_perV   = fitfunc_Coulomb_wF
        Params_perV = np.ones((2, Nconf))
    else:
        FuncName_perV, Params_perV = input_params(iFname_perV)
        if (FuncName_perV is Params_perV is None):
            quit()
        if (Nconf != len(Params_perV[0, :])):
            print("\nERROR: #.conf is differ between Wave[i] and per.Pot[i], exit.\n"); quit()
        
        Func_perV = set_fitfunc_from_fname(FuncName_perV)
        if (Func_perV is None):
            quit()
    
    if (SFP_flg):
        factor_perV = tmpRho0
        if (iFname_perV == 'Coulomb'):
            factor_perV = tmpRho0 * Np / tmpNA
        
        tmpParams_perV = Params_perV
        tmpFunc_perV   = Func_perV
        
        Func_perV   = lambda r,Ag1,Ag2,Ag3,Ag4: factor_perV * solve_Fint_1D(r,Ag1,Ag2,Ag3,Ag4)
        Params_perV = np.array([[None]*Nconf]*4)
        
        for i in range(Nconf):
            Params_perV[0, i] = integrand_1D_PotOpt_DensWS
            Params_perV[1, i] = tmpFunc_perV
            Params_perV[2, i] = tmpParams_perV[:, i]
            Params_perV[3, i] = Params_dens
    
    #print factor_perV; print Func_perV; print Params_perV; quit()
    
### Calculation & Output
    if (iFname_orgV is not None):
        Energy_0 = np.empty(Nconf)
        print("#")
        for i in range(Nconf):
            Energy_0[i] = solve_sch_0th_pert(Range, Coeff[:, i], mass, Func_orgV, Params_orgV[:, i])
            print("# conf = %d: End" % i)
        
        Mean, Err = make_mean_err(Energy_0)
        print("#\n# E_0 = %1.16f +/- %1.16f" % (Mean, Err))
    
    Energy_1 = np.empty(Nconf)
    print("#")
    for i in range(Nconf):
        Energy_1[i] = solve_sch_1st_pert(Range, Coeff[:, i], Func_perV, Params_perV[:, i])
        print("# conf = %d: End" % i)
    
    Mean, Err = make_mean_err(Energy_1)
    print("#\n# E_1 = %1.16f +/- %1.16f" % (Mean, Err))
    
    print("#\n# JOB END: elapsed_time = %lf [sec]" % (time.time() - JOB_start))
