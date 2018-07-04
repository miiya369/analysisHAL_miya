#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
hbar_c = 197.327053

if_V_C = None
if_V_T = None
ifengy = None
ofname = '/dev/null'
mass1  = 1000.0
mass2  = 1000.0
Emin   = 0.0
Edel   = 1.0
Emax   = 100.0
Rmax   = 10.0
Nproc  = 1

def fproc_wrapper(args):
    from sch_diffeq.solve_diffeq import solve_sch_diff
    a_iconf  = args[0]
    a_iphi   = args[1]
    a_ffuncs = args[2]
    a_params = args[3]
    a_mass   = args[4]
    a_Edata  = args[5]
    a_l      = args[6]
    a_Rini   = args[7]
    a_Rmax   = args[8]
    results = solve_sch_diff(a_iphi, a_ffuncs, a_params, a_mass, a_Edata, a_l, a_Rini, a_Rmax)
    print("# Calculate T-matrix... %3d end" % a_iconf)
    return results
### =========================== Main =========================== ###

def main():
    from multiprocessing import Pool
    
    from common.statistics             import make_mean_err
    from fitting.io_params             import input_params
    from fitting.fitfunc_type          import set_fitfunc_from_fname
    from sch_diffeq.solve_diffeq       import solve_sch_diff
    from Tmatrix.io_Tmatrix            import output_Tmatrix
    from Tmatrix.convert_mat           import convert_TtoS
    from Tmatrix.calc_phase_shift      import calc_phase_Sii, within_one
    
### Input the potential func_type and parameters ###
    func_name_V_C, params_V_C = input_params(if_V_C)
    if (func_name_V_C is params_V_C is None):
        return -1
    func_name_V_T, params_V_T = input_params(if_V_T)
    if (func_name_V_T is params_V_T is None):
        return -1
    mass  = np.array([[mass1, mass2], [mass1, mass2]])
    Nconf = len(params_V_C[:,0])
    
### Set the potential func_type ###
    fit_func_V_C = set_fitfunc_from_fname(func_name_V_C)
    fit_func_V_T = set_fitfunc_from_fname(func_name_V_T)
    
    fit_func = np.array([[lambda r,i: fit_func_V_C(r,*params_V_C[i,:]),
                          lambda r,i: fit_func_V_T(r,*params_V_T[i,:]) * np.sqrt(8)],
                         [lambda r,i: fit_func_V_T(r,*params_V_T[i,:]) * np.sqrt(8),
                          lambda r,i: fit_func_V_C(r,*params_V_C[i,:]) - fit_func_V_T(r,*params_V_T[i,:]) * 2.0]])
    params = np.array([[[[iconf] for j in range(2)] for i in range(2)] for iconf in range(Nconf)])
    
### Set the energies to calculate ###
    if (ifengy is None):
        N_E   = int((Emax - Emin) / Edel)
        Edata = np.array([Emin + i*Edel for i in range(N_E)])
    else:
        Edata = np.array([float(line.split()[0].strip()) for line in open(ifengy, 'r')])
        N_E   = len(Edata)
    
### Set the initial wave functions ###
    rini = 1e-6
    iphi = np.array([[rini**1, rini**3, 1.0, 3*rini**2],  # S-wave
                     [rini**1, rini**3, 0.0, 3*rini**2]]) # D-wave
    
### T-matrix calculation ###
    
    l_in  = np.array([0, 2])
    print("#\n# Calculate T-matrix...")
    if (Nproc == 1):
        Tmat = np.array([fproc_wrapper((iconf, iphi, fit_func, params[iconf,:], 
                                        mass, Edata, l_in, rini, Rmax)) for iconf in range(Nconf)])
    else:
        args_procs = [(iconf, iphi, fit_func, params[iconf,:], mass, Edata, l_in, rini, Rmax) for iconf in range(Nconf)]
        with Pool(Nproc) as proc:
            Tmat = np.array(proc.map(fproc_wrapper, args_procs))
    print("# Calculate T-matrix... all end\n#")
    
### Output T-matrix ###
    output_Tmatrix(ofname, Edata, Tmat)
    
### Output the Phase shift & Calculation of the Scattering length ###
    Smat = np.array([[convert_TtoS(Tmat[iconf,iE,:,:]) for iE in range(N_E)] for iconf in range(Nconf)])
    
    Edata[Edata==0.0] = 1e-10 # To avoid zero-div
    phase  = np.array([[calc_phase_Sii(Smat[iconf,:,i,i], shift_disc = 50)     for i  in range(2)  ] for iconf in range(Nconf)])
    MixAng = np.array([[np.arccos(within_one(phase[iconf,0,1,iE])) *90.0/np.pi for iE in range(N_E)] for iconf in range(Nconf)])
    Edata[Edata==1e-10] = 0.0
    
    print("#\n# E, phs(S), phs_e(S), phs(D), phs_e(D), mix_ang, mix_ang_e")
    for iE in range(N_E):
        print("%lf %e %e %e %e %e %e" % (Edata[iE],
                                         *make_mean_err(phase [:,0,0,iE]),
                                         *make_mean_err(phase [:,1,0,iE]),
                                         *make_mean_err(MixAng[:,    iE])))
    
    # Scattering length for [fm] := T_ii(p)/p (p->0)
    ScatLength = np.array([solve_sch_diff(iphi, fit_func, params[iconf,:],
                                          mass, np.array([1e-10]), l_in, rini, Rmax)[0,:,:]
                           for iconf in range(Nconf)]) / np.sqrt(2.0 * mass * 1e-10) * hbar_c
    print("#\n# Scattering length for S-wave[fm] = %lf(%lf) + %lf(%lf) i" % (*make_mean_err(ScatLength[:,0,0].real),
                                                                             *make_mean_err(ScatLength[:,0,0].imag)))
    print("#\n# Scattering length for D-wave[fm] = %lf(%lf) + %lf(%lf) i" % (*make_mean_err(ScatLength[:,1,1].real),
                                                                             *make_mean_err(ScatLength[:,1,1].imag)))
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifparams (V_C)] [ifparams (V_T)] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --ofile [Output file name (For T-matrix)   ] Default =", ofname)
    print("      --mass1 [Mass of 1st baryon                ] Default =", mass1)
    print("      --mass2 [Mass of 2nd baryon                ] Default =", mass2)
    print("      --Efile [Input file name (For Energy cood.)] Default =", ifengy)
    print("      --Emin  [Minimum energy  for output        ] Default =", Emin)
    print("      --Edel  [energy division for output        ] Default =", Edel)
    print("      --Emax  [Maximum energy  for output        ] Default =", Emax)
    print("      --Rmax  [Maximum  range  for output        ] Default =", Rmax)
    print("      --Nproc [#.process                         ] Default =", Nproc)
    exit("\nNote: The potentials, masses, energies and range are should be physical unit (MeV, fm).")

def check_args():
    print("# === Check Arguments ===")
    print("# if VC =", if_V_C)
    print("# if VT =", if_V_T)
    print("# fengy =", ifengy)
    print("# ofile =", ofname)
    print("# mass1 =", mass1)
    print("# mass2 =", mass2)
    print("# E min =", Emin)
    print("# E del =", Edel)
    print("# E max =", Emax)
    print("# R max =", Rmax)
    print("# Nproc =", Nproc)
    print("# =======================")

def set_args(ARGC, ARGV):
    global if_V_C, if_V_T, ifengy, ofname, mass1, mass2, Emin, Edel, Emax, Rmax, Nproc
    
    if (ARGC < 3 or ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    if_V_C = ARGV[1].strip()
    if_V_T = ARGV[2].strip()
    
    for i in range(3, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass1'):
                mass1 = float(ARGV[i+1])
            elif (ARGV[i] == '--mass2'):
                mass2 = float(ARGV[i+1])
            elif (ARGV[i] == '--ofile'):
                ofname = ARGV[i+1].strip()
            elif (ARGV[i] == '--Efile'):
                ifengy = ARGV[i+1].strip()
                Emin   = 'From the file'
                Edel   = 'From the file'
                Emax   = 'From the file'
            elif (ARGV[i] == '--Emin'):
                if (ifengy is not None):
                    print("\nWARNING: '--Emin' option is ignored " +
                          "when '--Efile' option was specified.\n")
                else:
                    Emin = float(ARGV[i+1])
            elif (ARGV[i] == '--Edel'):
                if (ifengy is not None):
                    print("\nWARNING: '--Edel' option is ignored " +
                          "when '--Efile' option was specified.\n")
                else:
                    Edel = float(ARGV[i+1])
            elif (ARGV[i] == '--Emax'):
                if (ifengy is not None):
                    print("\nWARNING: '--Emax' option is ignored " +
                          "when '--Efile' option was specified.\n")
                else:
                    Emax = float(ARGV[i+1])
            elif (ARGV[i] == '--Rmax'):
                Rmax = float(ARGV[i+1])
            elif (ARGV[i] == '--Nproc'):
                Nproc = int(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
    check_args()
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
