#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname1   = None
ifname2   = None
ofname    = "/dev/null"
func_name = "2SG"
min_range = 0.0
max_range = 3.0
Nrestart  = 1
fparam    = None
iparams   = np.empty(0)
Upper_bnd = np.empty(0)
Lower_bnd = np.empty(0)

Max_iter  = 100000
Nparam    = 0

def chisq_i(p, f, x, y_n, y_d, e):
    return (f(x,*p)*y_d - y_n) / e
### =========================== Main =========================== ###

def main ():
    global iparams, Upper_bnd, Lower_bnd
    
    from scipy.optimize       import leastsq
    from scipy.optimize       import least_squares
    
    from common.io_data_bin   import input_bin_data
    from common.statistics    import make_mean_err
    from fitting.io_params    import output_params
    from fitting.fitfunc_type import set_fitfunc_from_fname
    from fitting.print_gnuplt import print_gnuplt
    
### Input data ###
    yData1_tmp, xData1_tmp, eData1_tmp = input_bin_data(ifname1)
    if (xData1_tmp is yData1_tmp is eData1_tmp is None):
        return -1;
    yData2_tmp, xData2_tmp, eData2_tmp = input_bin_data(ifname2)
    if (xData2_tmp is yData2_tmp is eData2_tmp is None):
        return -1;
    
    Nconf      = len(yData1_tmp[:,0])
    Ndata_tmp  = len(yData1_tmp[0,:])
    Nparam     = len(iparams)
    
    if (Nconf     != len(yData2_tmp[:,0])):
        print("\nERROR: #.conf in the two files are differ, exit."); return -1
    if (Ndata_tmp != len(yData2_tmp[0,:])):
        print("\nERROR: #.data in the two files are differ, exit."); return -1
    
    Ndata = 0
    yData1 = np.empty((Nconf,0))
    yData2 = np.empty((Nconf,0))
    xData  = np.empty(0)
    eData1 = np.empty(0)
    eData2 = np.empty(0)
    for ir in range(Ndata_tmp):
        if (xData1_tmp[ir] != xData2_tmp[ir]):
            print("\nERROR: x-data in the two files are differ, exit."); return -1
        if (min_range <= xData1_tmp[ir] <= max_range):
            yData1 = np.append(yData1, yData1_tmp[:,ir].reshape((Nconf,1)), axis=1)
            yData2 = np.append(yData2, yData2_tmp[:,ir].reshape((Nconf,1)), axis=1)
            xData  = np.append(xData , xData1_tmp[  ir])
            eData1 = np.append(eData1, eData1_tmp[  ir])
            eData2 = np.append(eData2, eData2_tmp[  ir])
            Ndata += 1
    print("#"); print_gnuplt(func_name, iparams); print("#")
    del yData1_tmp
    del xData1_tmp
    del eData1_tmp
    del yData2_tmp
    del xData2_tmp
    del eData2_tmp
    
### Take fit ###
    fit_func = set_fitfunc_from_fname(func_name)
    if (np.all(Lower_bnd == -np.inf) and np.all(Upper_bnd == np.inf)):
        LU_bounds = None
    else:
        LU_bounds = (Lower_bnd[:], Upper_bnd[:])
    
    for ires in range(Nrestart):
        if (LU_bounds is None):
            res_params = np.array([leastsq(chisq_i, iparams, 
                                           args=(fit_func, xData, yData1[iconf,:], yData2[iconf,:], eData1),
                                           maxfev=Max_iter)[0] for iconf in range(Nconf)])
        else:
            res_params = np.array([least_squares(chisq_i, iparams, 
                                                 args=(fit_func, xData, yData1[iconf,:], yData2[iconf,:], eData1),
                                                 bounds=LU_bounds, max_nfev=Max_iter).x for iconf in range(Nconf)])
        
        chisq_dof = np.array([np.sum(chisq_i(res_params[iconf,:], fit_func, xData, 
                                             yData1[iconf,:], yData2[iconf,:], eData1)**2) /
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
###### Functions for arguments ######
def usage(ARGV0):
    print("usage  : python %s [ifname1 (data for numelator)] [ifname2 (data for denominator)] {options}\n" % 
          os.path.basename(ARGV0))
    print("brief  : Minimize sum[{(y_d * f(x,*p) - y_n) / (e_n)}**2],")
    print("       : instead of sum[{(f(x,*p) - y) / e}**2], where y = y_n/y_d.\n")
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
    print("# ifile 1   =", ifname1)
    print("# ifile 2   =", ifname2)
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
    
    global ifname1, ifname2, ofname, func_name, min_range, max_range, Nrestart
    global iparams, fparam, Upper_bnd, Lower_bnd, Nparam
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname1 = ARGV[1].strip()
    ifname2 = ARGV[2].strip()
    
    for i in range(3, ARGC):
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
    if (argc < 3):
        usage(argv[0])

    set_args(argc, argv)

    t_start = time.time()
    if (main() != 0):
        exit("ERROR EXIT.")
    print("#\n# Elapsed time [s] = %d" % (time.time() - t_start))
