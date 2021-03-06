#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname    = None
ofname    = "/dev/null"
func_name = "3G"
min_range = 0.001
max_range = 2.0
Nrestart  = 1
fparam    = None
iparams   = np.empty(0)
Upper_bnd = np.empty(0)
Lower_bnd = np.empty(0)

Max_iter  = 100000
Nparam    = 0
### =========================== Main =========================== ###

def main ():
    global iparams, Upper_bnd, Lower_bnd
    
    from scipy.optimize       import curve_fit
    
    from common.io_data_bin   import input_bin_data
    from common.statistics    import make_mean_err
    from fitting.io_params    import output_params
    from fitting.fitfunc_type import set_fitfunc_from_fname
    from fitting.print_gnuplt import print_gnuplt
    
### Input data ###
    yData_tmp, xData_tmp, eData_tmp = input_bin_data(ifname)
    if (xData_tmp is yData_tmp is eData_tmp is None):
        return -1;
    
    Nconf      = len(yData_tmp[:,0])
    Ndata_tmp  = len(yData_tmp[0,:])
    Nparam     = len(iparams)
    
    Ndata = 0
    yData = np.empty((Nconf,0))
    xData = np.empty(0)
    eData = np.empty(0)
    for ir in range(Ndata_tmp):
        if (min_range <= xData_tmp[ir] <= max_range):
            yData = np.append(yData, yData_tmp[:,ir].reshape((Nconf,1)), axis=1)
            xData = np.append(xData, xData_tmp[  ir])
            eData = np.append(eData, eData_tmp[  ir])
            Ndata += 1
    print("#"); print_gnuplt(func_name, iparams); print("#")
    del yData_tmp
    del xData_tmp
    del eData_tmp
    
### Take fit ###
    fit_func = set_fitfunc_from_fname(func_name)
    if (np.all(Lower_bnd == -np.inf) and np.all(Upper_bnd == np.inf)):
        LU_bounds = None
    else:
        LU_bounds = (Lower_bnd[:], Upper_bnd[:])
    
    for ires in range(Nrestart):
        if (LU_bounds is None):
            res_params = np.array([curve_fit(fit_func, xData, yData[iconf,:], 
                                             sigma=eData, p0=iparams, maxfev=Max_iter)[0]
                                   for iconf in range(Nconf)])
        else:
            res_params = np.array([curve_fit(fit_func, xData, yData[iconf,:], 
                                             sigma=eData, p0=iparams, bounds=LU_bounds, max_nfev=Max_iter)[0]
                                   for iconf in range(Nconf)])
        
        chisq_dof = np.array([np.sum(((fit_func(xData, *res_params[iconf,:]) - yData[iconf,:]) / eData)**2) /
                              (Ndata - Nparam) for iconf in range(Nconf)])
        
        print("# === Fitting Results (%03d) ===" % (ires+1))
        print("# chisq/dof = %lf +/- %lf" % (*make_mean_err(chisq_dof),))
        iparams = np.array([make_mean_err(res_params[:,iparam])[0] for iparam in range(Nparam)])
        print("# Parameters:", "%e "*len(iparams) % (*iparams,))
    
### Print Results ###
    print("#\n# === Fitting Results (Fin) ===")
    mean, err = make_mean_err(chisq_dof)
    print("# chisq/dof = %15lf +/- %15lf (%15.6f %%)" % (mean, err, abs(err/mean) * 100))
    
    tmp_params = np.array([make_mean_err(res_params[:,iparam]) for iparam in range(Nparam)])
    for iparam in range(Nparam):
        print("# param[%2d] = %15e +/- %15e (%15.6f %%)" % (iparam, tmp_params[iparam,0], tmp_params[iparam,1],
                                                            abs(tmp_params[iparam,1] / tmp_params[iparam,0]) * 100))
    
    print("#"); print_gnuplt(func_name, tmp_params[:,0]); print("#")
    output_params(ofname, func_name, res_params)
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile (miyamoto-format binary file)] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --ofile    [Output (parameters) file name               ] Default =", ofname)
    print("      --fit_func [Fit function name                           ] Default =", func_name)
    print("      --min_r    [Minimum range for fitting                   ] Default =", min_range)
    print("      --max_r    [Maximum range for fitting                   ] Default =", max_range)
    print("      --Nres     [#.restart of fitting                        ] Default =", Nrestart)
    print("      --fparam   [Initial parameters file name                ] Default =", fparam)
    print("      --params   [Initial parameters      (separated by space)] Default =  1.0 for all parameters")
    print("      --lower    [Lower bounds for params (separated by space)] Default = -inf for all parameters")
    print("      --upper    [Upper bounds for params (separated by space)] Default =  inf for all parameters")
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     =", ifname)
    print("# ofile     =", ofname)
    print("# fit func  =", func_name)
    print("# min range =", min_range)
    print("# max range =", max_range)
    print("# N.restart =", Nrestart)
    print("# N.param   =", Nparam)
    print("# ini.param =", iparams)
    print("# lower bnd =", Lower_bnd)
    print("# upper bnd =", Upper_bnd)
    print("# =======================")

def set_args(ARGC, ARGV):
    from fitting.io_params    import input_init_params
    from fitting.fitfunc_type import get_Nparam_from_fname
    
    global ifname, ofname, func_name, min_range, max_range, Nrestart
    global iparams, fparam, Upper_bnd, Lower_bnd, Nparam
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--ofile'):
                ofname    = ARGV[i+1].strip()
            elif (ARGV[i] == '--fit_func'):
                func_name  = ARGV[i+1].strip()
            elif (ARGV[i] == '--min_r'):
                min_range = float(ARGV[i+1])
            elif (ARGV[i] == '--max_r'):
                max_range = float(ARGV[i+1])
            elif (ARGV[i] == '--Nres'):
                Nrestart  = int(ARGV[i+1])
            elif (ARGV[i] == '--params'):
                iparams   = np.empty(0)
                for iparam in range(i+1, ARGC):
                    try:
                        iparams = np.append(iparams, float(ARGV[iparam]))
                    except:
                        break
            elif (ARGV[i] == '--upper'):
                Upper_bnd = np.empty(0)
                for iparam in range(i+1, ARGC):
                    try:
                        Upper_bnd = np.append(Upper_bnd, float(ARGV[iparam]))
                    except:
                        break
            elif (ARGV[i] == '--lower'):
                Lower_bnd = np.empty(0)
                for iparam in range(i+1, ARGC):
                    try:
                        Lower_bnd = np.append(Lower_bnd, float(ARGV[iparam]))
                    except:
                        break
            elif (ARGV[i] == '--fparam'):
                fparam    = ARGV[i+1].strip()
            else:
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
    if (get_Nparam_from_fname(func_name) is None):
        exit(1)
    else:
        Nparam = get_Nparam_from_fname(func_name)
    
    if (fparam is not None):
        iparams, Lower_bnd, Upper_bnd = input_init_params(fparam)
        if (iparams is Upper_bnd is Lower_bnd is None):
            exit(1)
    else:
        if (len(iparams  ) == 0):
            for iparam in range(Nparam):
                iparams   = np.append(iparams, 1.0)
        if (len(Upper_bnd) == 0):
            for iparam in range(Nparam):
                Upper_bnd = np.append(Upper_bnd, np.inf)
        if (len(Lower_bnd) == 0):
            for iparam in range(Nparam):
                Lower_bnd = np.append(Lower_bnd, -np.inf)
    
    check_args()
    
    if (len(iparams) != Nparam or len(Upper_bnd) != Nparam or len(Lower_bnd) != Nparam):
        exit("\n#.initial parameter is differ from the one of fit function, exit.\n")
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
