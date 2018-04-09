#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np

from MiscFuncs.Misc           import *
from MiscFuncs.IO_BinData     import *
from MiscFuncs.DataStatistics import make_mean_err
from Fitting.SimulFit         import *
from Fitting.IO_Params        import *
from Fitting.PrintGnuForm     import *

iFname1   = None
iFname2   = None
oFname    = "/dev/null"
Fparam    = None
fitF_name = "3G_Simul"
min_range = 0
max_range = 25
in_Params = np.empty(0)
Nrestart  = 1
Upper_b   = np.empty(0)
Lower_b   = np.empty(0)

Fparam_flg = False

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [input bin-data 1] [input bin-data 2] {--option}\n" % ARGV[0])
    print("option:")
    print("      --ofile    [Output (parameters) file name          ] Default ="),; print oFname
    print("      --fit_func [Fit function name                      ] Default ="),; print fitF_name
    print("      --min_r    [Minimum range for fitting              ] Default ="),; print min_range
    print("      --max_r    [Maximum range for fitting              ] Default ="),; print max_range
    print("      --Nres     [#.restart of fitting                   ] Default ="),; print Nrestart
    print("      --fparam   [Initial parameters file name           ] Default ="),; print Fparam
    print("      --params   [Initial parameters      (separated by space)] Default =  1.0 for all parameters")
    print("      --lower    [Lower bounds for params (separated by space)] Default = -inf for all parameters")
    print("      --upper    [Upper bounds for params (separated by space)] Default =  inf for all parameters")
    print; quit()

def set_args(ARGC, ARGV):
    global iFname1, iFname2, oFname, fitF_name, min_range, max_range, Nrestart, in_Params, Fparam, Fparam_flg
    global Upper_b, Lower_b
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFname1 = ARGV[1].strip()
    iFname2 = ARGV[2].strip()
    
    for i in range(3, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--ofile'):
                oFname = ARGV[i+1].strip()
            elif (ARGV[i] == '--fit_func'):
                fitF_name = ARGV[i+1].strip()
            elif (ARGV[i] == '--min_r'):
                min_range = float(ARGV[i+1])
            elif (ARGV[i] == '--max_r'):
                max_range = float(ARGV[i+1])
            elif (ARGV[i] == '--Nres'):
                Nrestart = int(ARGV[i+1])
            elif (ARGV[i] == '--params'):
                Fparam_flg = False; in_Params = np.empty(0)
                for iparam in range(i+1, ARGC):
                    try:
                        in_Params = np.append(in_Params, float(ARGV[iparam]))
                    except:
                        break
            elif (ARGV[i] == '--upper'):
                Fparam_flg = False; Upper_b = np.empty(0)
                for iparam in range(i+1, ARGC):
                    try:
                        Upper_b = np.append(Upper_b, float(ARGV[iparam]))
                    except:
                        break
            elif (ARGV[i] == '--lower'):
                Fparam_flg = False; Lower_b = np.empty(0)
                for iparam in range(i+1, ARGC):
                    try:
                        Lower_b = np.append(Lower_b, float(ARGV[iparam]))
                    except:
                        break
            elif (ARGV[i] == '--fparam'):
                Fparam_flg = True
                Fparam     = ARGV[i+1].strip()
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i]); usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile 1   ="),; print iFname1
    print("# ifile 2   ="),; print iFname2
    print("# ofile     ="),; print oFname
    print("# fit func  ="),; print fitF_name
    print("# min range ="),; print min_range
    print("# max range ="),; print max_range
    print("# N.restart ="),; print Nrestart
    print("# ini.param ="),; print in_Params
    print("# lower bnd ="),; print Lower_b
    print("# upper bnd ="),; print Upper_b
    print("# =======================")

###### Main part
if __name__ == "__main__":
    from scipy.optimize import curve_fit
    
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv)
    
    if (Fparam_flg):
        in_Params, Lower_b, Upper_b = input_init_params(Fparam)
        if (in_Params is Upper_b is Lower_b is None):
            quit()
    else:
        if (len(in_Params) == 0):
            for iparam in range(get_Nparam_from_fname_simul(fitF_name)):
                in_Params = np.append(in_Params, 1.0)
        if (len(Upper_b) == 0):
            for iparam in range(get_Nparam_from_fname_simul(fitF_name)):
                Upper_b = np.append(Upper_b, np.inf)
        if (len(Lower_b) == 0):
            for iparam in range(get_Nparam_from_fname_simul(fitF_name)):
                Lower_b = np.append(Lower_b, -np.inf)
    
    check_args()
    
    if (len(in_Params) != get_Nparam_from_fname_simul(fitF_name) or
        len(Upper_b)   != get_Nparam_from_fname_simul(fitF_name) or
        len(Lower_b)   != get_Nparam_from_fname_simul(fitF_name)):
        print("\n#.initial parameter is differ from the one of fit function, exit.\n"); quit()
    
### Input data ###
    yData_tmp1, xData_tmp1, eData_tmp1 = input_bin_data(iFname1)
    if (xData_tmp1 is yData_tmp1 is eData_tmp1 is None):
        quit()
    yData_tmp2, xData_tmp2, eData_tmp2 = input_bin_data(iFname2)
    if (xData_tmp2 is yData_tmp2 is eData_tmp2 is None):
        quit()
    
    if (len(yData_tmp1[0, :]) != len(yData_tmp2[0, :])):
        print("\nERROR: Different #.conf, exit.\n"); quit()
    if (len(yData_tmp1[:, 0]) != len(yData_tmp2[:, 0])):
        print("\nERROR: Different #.data, exit.\n"); quit()
    for idata in range(len(yData_tmp1[:, 0])):
        if (xData_tmp1[idata] != xData_tmp2[idata]):
            print("\nERROR: Different x-data, exit.\n"); quit()
    
    Nconf      = len(yData_tmp1[0, :])
    Ndata_tmp  = len(yData_tmp1[:, 0])
    Nparam     = len(in_Params)
    
    Ndata = 0
    yData = np.empty((0,Nconf))
    xData = np.empty(0)
    eData = np.empty(0)
    for ir in range(Ndata_tmp):
        if (min_range <= xData_tmp1[ir] <= max_range):
            yData = np.append(yData, np.empty((1,Nconf)), axis=0)
            for iconf in range(Nconf):
                yData[Ndata, iconf] = yData_tmp1[ir, iconf]
            xData = np.append(xData, xData_tmp1[ir])
            eData = np.append(eData, eData_tmp1[ir])
            
            Ndata += 1
            
            yData = np.append(yData, np.empty((1,Nconf)), axis=0)
            for iconf in range(Nconf):
                yData[Ndata, iconf] = yData_tmp2[ir, iconf]
            xData = np.append(xData, xData_tmp2[ir] + 100)
            eData = np.append(eData, eData_tmp2[ir])
            
            Ndata += 1
    
    print("#")
    tmp_fitF_name = convert_simul_fitparam(in_Params, fitF_name, 0)[0]
    print_gnu(tmp_fitF_name, convert_simul_fitparam(in_Params, fitF_name, 0)[1])
    print_gnu(tmp_fitF_name, convert_simul_fitparam(in_Params, fitF_name, 1)[1])
    print("#")
    del yData_tmp1
    del xData_tmp1
    del eData_tmp1
    del yData_tmp2
    del xData_tmp2
    del eData_tmp2
    
### Take fit ###
    ResParams = np.empty((Nparam, Nconf))
    Chisq_dof = np.empty(Nconf)
    fit_func  = set_fitfunc_from_fname_simul(fitF_name)
    LU_bounds = (Lower_b[:], Upper_b[:])
    
    for ires in range(Nrestart):
        for iconf in range(Nconf):
            tmpParams, cov   = curve_fit(fit_func, xData, yData[:, iconf], 
                                         sigma=eData, p0=in_Params, bounds=LU_bounds)
                                         #, maxfev=Max_iter)
                                         #, max_nfev=Max_iter)
            
            Chisq_dof[iconf] = np.sum(((fit_func(xData, *tmpParams) - yData[:, iconf]) 
                                       / eData)**2) / (Ndata - Nparam)
            
            for iparam in range(Nparam):
                ResParams[iparam, iconf] = tmpParams[iparam]
        
        print("# === Fitting Results (%d) ===" % (ires+1))
        tmp_mean, tmp_err = make_mean_err(Chisq_dof)
        print("# Chisq/dof = %lf +/- %lf" % (tmp_mean, tmp_err))
        print("# Parameters:"),
        for iparam in range(Nparam):
            tmp_mean, tmp_err = make_mean_err(ResParams[iparam, :])
            in_Params[iparam] = tmp_mean
            print("%e" % tmp_mean),
        print
    
### Print Results ###
    print("#\n# === Fitting Results (Final) ===")
    tmp_mean, tmp_err = make_mean_err(Chisq_dof)
    print("# Chisq/dof = %lf +/- %lf" % (tmp_mean, tmp_err))
    tmpParams = np.empty(Nparam)
    for iparam in range(Nparam):
        tmpParams[iparam], tmp_err = make_mean_err(ResParams[iparam, :])
        print("# param[%2d] = %e +/- %e" % (iparam, tmpParams[iparam], tmp_err))
    
    print("#")
    tmp_fitF_name = convert_simul_fitparam(tmpParams, fitF_name, 0)[0]
    print_gnu(tmp_fitF_name, convert_simul_fitparam(tmpParams, fitF_name, 0)[1])
    print_gnu(tmp_fitF_name, convert_simul_fitparam(tmpParams, fitF_name, 1)[1])
    print("#")
    
    out_params(oFname, fitF_name, ResParams)
