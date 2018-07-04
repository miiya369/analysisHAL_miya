#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc                    import *
from MiscFuncs.DataStatistics          import make_mean_err
from Fitting.IO_Params                 import input_params, out_params
from FoldingPotential.DensityFunc      import *
from FoldingPotential.Integrand        import *
from FoldingPotential.SolveIntegral    import *
from FoldingPotential.PotentialFunc    import *
from FoldingPotential.ConvertFitParams import *
from Fitting.FitFunctionType           import *
from Fitting.FitFunctionForm           import fitfunc_Coulomb, fitfunc_Coulomb_wF

iFname    = None
N_p       = 82
N_A       = 208
R_A       = 6.80
dif_A     = 0.515
Nplot     = 100
max_r     = 12
dens_type = "Woods_Saxon"
with_Coul = False
oFname    = None

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [Fit Param] {--option}\n" % ARGV[0])
    print("option:")
    print("      --Np          [#.proton (For Coulomb)] Default ="),; print N_p
    print("      --A           [#.nucleus             ] Default ="),; print N_A
    print("      --RA          [Radias of nuclei  (fm)] Default ="),; print R_A
    print("      --diffuse     [Diffuseness       (fm)] Default ="),; print dif_A
    print("      --Nplot       [#.plot for output     ] Default ="),; print Nplot
    print("      --max_r       [max r  for output (fm)] Default ="),; print max_r
    print("      --dtype       [Density type          ] Default ="),; print dens_type
    print("      --conv_param  [Output file name] <- For ONLY Gaussian-type density"),; print
    print("      --w_Coulomb    Calculate with Coulomb potential"),; print
    print; quit()

def set_args(ARGC, ARGV):
    global iFname, N_p, N_A, R_A, dif_A, Nplot, max_r, dens_type, oFname, with_Coul
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--A'):
                N_A = int(ARGV[i+1])
            elif (ARGV[i] == '--RA'):
                R_A = float(ARGV[i+1])
            elif (ARGV[i] == '--Np'):
                N_p = int(ARGV[i+1])
            elif (ARGV[i] == '--diffuse'):
                dif_A = float(ARGV[i+1])
            elif (ARGV[i] == '--Nplot'):
                Nplot = int(ARGV[i+1])
            elif (ARGV[i] == '--max_r'):
                max_r = float(ARGV[i+1])
            elif (ARGV[i] == '--dtype'):
                dens_type = ARGV[i+1].strip()
            elif (ARGV[i] == '--conv_param'):
                oFname = ARGV[i+1].strip()
            elif (ARGV[i] == '--w_Coulomb'):
                with_Coul = True
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     ="),; print iFname
    print("# N_p       ="),; print N_p
    print("# A         ="),; print N_A
    print("# R_A       ="),; print R_A
    print("# diffuse   ="),; print dif_A
    print("# N.plot    ="),; print Nplot
    print("# max range ="),; print max_r
    print("# dens.type ="),; print dens_type
    print("# ofile     ="),; print oFname
    print("# w/Coulomb ="),; print with_Coul
    print("# =======================")

###### Main part
if __name__ == "__main__":
    import numpy as np
    import re
    
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv)
    check_args()
    
### Input parameters ###
    FuncName_pot, Params_pot = input_params(iFname)
    if (FuncName_pot is Params_pot is None):
        quit()
    
    Nconf = len(Params_pot[0, :])
    
    Func_pot = set_fitfunc_from_fname(FuncName_pot)
    if (Func_pot is None):
        quit()
    
### Construct the Folding potential & Output ###
    if   (dens_type == 'Woods_Saxon'):
        Integrand_1D = integrand_1D_PotOpt_DensWS
        Func_dens    = dens_woods_saxon
        Params_dens  = np.array((R_A, dif_A))
    elif (dens_type == 'Gaussian'):
        Integrand_1D = integrand_1D_PotOpt_DensGauss
        Func_dens    = dens_gaussian
        Params_dens  = np.array((R_A,))
    else:
        print("\nERROR: Invalid density type: %s, exit.\n" % dens_type)
        quit()
    
    rho_0 = calc_rho0(N_A, Func_dens, Params_dens)
    if (rho_0 is None):
        quit()
    print("#\n# rho_0 [fm^{-3}] ="),; print rho_0; print("#")
    
    Fpot_V = lambda r,Ag1,Ag2,Ag3,Ag4: solve_Fint_1D(r,Ag1,Ag2,Ag3,Ag4) * rho_0
    Fpot_C = lambda r,Ag1,Ag2,Ag3,Ag4: solve_Fint_1D(r,Ag1,Ag2,Ag3,Ag4) * rho_0 * N_p / N_A
    
    if (oFname is None):
        Results = np.empty(Nconf);
        factor  = max_r / float(Nplot)
        
        for iplot in range(1, Nplot+1):
            r_fpot = iplot * factor;
            
            for iconf in range(Nconf):
                Results[iconf] = Fpot_V(r_fpot, Integrand_1D, Func_pot,
                                        Params_pot[:, iconf], Params_dens)
            
            if (with_Coul):
                if (dens_type == 'Woods_Saxon'):
                    Results_Coul_woF = fpot_Coulomb_dens_WS(r_fpot, rho_0, R_A, dif_A, N_p, N_A)
                    Results_Coul_w_F = Fpot_C(r_fpot, Integrand_1D, fitfunc_Coulomb_wF,
                                              np.array((1.0, 1.0)), Params_dens)
                else:
                    Results_Coul_woF = Fpot_C(r_fpot, Integrand_1D, fitfunc_Coulomb,
                                              np.array((1.0, 1.0)), Params_dens)
                    Results_Coul_w_F = Fpot_C(r_fpot, Integrand_1D, fitfunc_Coulomb_wF,
                                              np.array((1.0, 1.0)), Params_dens)
                
                #print("%lf %1.16lf %1.16lf" % (r_fpot, Results_Coul_woF, Results_Coul_w_F))
                #Results += Results_Coul_woF
                Results += Results_Coul_w_F
            
            Mean, Err = make_mean_err(Results)
            
            print("%lf %1.16e %1.16e" % (r_fpot, Mean, Err))
            
            #print("%lf" % r_fpot), # For all configurations
            #for iconf in range(Nconf):
            #    print(" %e" % Results[iconf]),
            #print
    
    else:
        if (dens_type != 'Gaussian'):
            print("\nERROR: The density-type '%s' cannot be converted, exit.\n" % dens_type)
            quit()
        if (re.match('^[1-9]G$', FuncName_pot) is None):
            print("\nERROR: Only the Multi-Gaussian-type potentials can be converted, exit.\n")
            quit()
        
        out_params(oFname, FuncName_pot, convert_params_gauss(Params_pot, N_A, R_A))
