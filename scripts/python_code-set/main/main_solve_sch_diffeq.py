#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np
import time

from common.statistics         import make_mean_err
from fitting.io_params         import input_params
from fitting.fitfunc_type      import set_fitfunc_from_fname
from sch_diffeq.io_param_files import input_NxN_params
from sch_diffeq.solve_diffeq   import solve_sch_diff
from Tmatrix.io_Tmatrix        import output_Tmatrix
from Tmatrix.convert_mat       import convert_TtoS
from Tmatrix.calc_phase_shift  import calc_phase_Nch1

### ================== Global Parameters Init. ================= ###
hbar_c = 197.327053

ifname = None
ifengy = None
ofname = '/dev/null'
mass1  = 1000.0
mass2  = 1000.0
Emin   = 0.0
Edel   = 1.0
Emax   = 100.0
Rmax   = 10.0
Nch    = 1
### =========================== Main =========================== ###
def main():
    t_start = time.time()
    
### Input the potential func_type and parameters ###
    if (Nch == 1):
        TmpFuncName, TmpParams = input_params(ifname)
        if (TmpFuncName is TmpParams is None):
            return -1
        FuncName = np.array([[TmpFuncName ]])
        Params   = np.array([[TmpParams   ]])
        mass     = np.array([[mass1, mass2]])
    else:
        FuncName, Params, mass = input_NxN_params(ifname)
        if (FuncName is Params is mass is None):
            return -1
        if (Nch != len(Params[:,0,0,0])):
            print("\nERROR: Different #.ch in the NxN parameters file " +
                  "(%d != %d), exit." % (Nch, len(Params[:,0,0,0]))); return -1
    
### Set the potential func_type ###
    fit_funcs = np.array([[set_fitfunc_from_fname(FuncName[ich][jch])
                           for jch in range(Nch)] for ich in range(Nch)])
    
### Set the energies to calculate ###
    if (ifengy is None):
        N_E   = int((Emax - Emin) / Edel)
        Edata = np.array([Emin + i*Edel for i in range(N_E)])
    else:
        Edata = np.array([float(line.split()[0].strip()) for line in open(ifengy, 'r')])
        N_E   = len(Edata)
    
### T-matrix calculation ###
    Nconf  = len(Params[0,0,:,0])
    Nparam = len(Params[0,0,0,:])
    Tmat   = np.empty((Nconf, N_E, Nch, Nch), dtype=complex)
    print("#\n# Calculate T-matrix...")
    for iconf in range(Nconf):
        tmpParams = np.array([[[Params[ich,jch,iconf,iparam]
                                for iparam in range(Nparam)] 
                               for jch in range(Nch)] for ich in range(Nch)])
        Tmat[iconf,:,:,:] = solve_sch_diff(fit_funcs, tmpParams, mass, Edata, Rmax)
        print("# Calculate T-matrix... end: iconf=%03d" % iconf)
    
### Output T-matrix ###
    output_Tmatrix(ofname, Edata, Tmat)
    if (Nch != 1):
        return 0
    
### Output the Phase shift & Calculation of the Scattering length (For #.ch == 1) ###
    convert_TtoS(Tmat)
    
    tmp_mu = (mass1*mass2)/(mass1+mass2)
    Edata[Edata==0.0] = 1e-10 # To avoid zero-div
    
    PhaseShift = np.array([[calc_phase_Nch1(Tmat[iconf,iE,0,0])[0]
                            for iE in range(N_E)] for iconf in range(Nconf)])
    
    ScatLength = np.array([[(np.tan (PhaseShift[iconf] * np.pi/180) /
                             np.sqrt(2.0 * tmp_mu * Edata[iE]) * hbar_c)
                            for iE in range(N_E)] for iconf in range(Nconf)])
    
    Edata[Edata==1e-10] = 0.0
    print("#\n# E, phs, phs_e, scatt.len, scatt.len_e")
    for iE in range(N_E):
        mean1, err1 = make_mean_err(PhaseShift[:,iE])
        mean2, err2 = make_mean_err(ScatLength[:,iE])
        print("%lf %e %e %e %e" % (Edata[iE], mean1, err1, mean2, err2))
    
    # Scattering length [fm] := T(p)/p (p->0)
    TmpTmat    = np.array([solve_sch_diff(fit_funcs, np.array([[Params[0,0,iconf,:]]]),
                                          mass, np.array([1e-10]), Rmax)[0,0,0]
                          for iconf in range(Nconf)])
    ScatLength = TmpTmat / np.sqrt(2.0 * tmp_mu * 1e-10) * hbar_c
    m_r, e_r   = make_mean_err(ScatLength.real)
    m_i, e_i   = make_mean_err(ScatLength.imag)
    print("#\n# Scattering length [fm] = %lf(%lf) + %lf(%lf) i" % (m_r, e_r, m_i, e_i))
    
    print("#\n# Elapsed time [s] = %d" % (time.time() - t_start))
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [params data or NxN summarized file] {options}" % os.path.basename(ARGV0))
    print("options:")
    print("      --ofile [Output file name (For T-matrix)   ] Default =", ofname)
    print("      --mass1 [Mass of 1st baryon (For #.ch=1)   ] Default =", mass1)
    print("      --mass2 [Mass of 2nd baryon (For #.ch=1)   ] Default =", mass2)
    print("      --Efile [Input file name (For Energy cood.)] Default =", ifengy)
    print("      --Emin  [Minimum energy  for output        ] Default =", Emin)
    print("      --Edel  [energy division for output        ] Default =", Edel)
    print("      --Emax  [Maximum energy  for output        ] Default =", Emax)
    print("      --Rmax  [Maximum  range  for output        ] Default =", Rmax)
    print("      --Nch   [#.channel                         ] Default =", Nch)
    exit("\nNote: The mass, energy and range are should be physical unit (MeV, fm).")

def check_args():
    print("# === Check Arguments ===")
    print("# ifile =", ifname)
    print("# fengy =", ifengy)
    print("# ofile =", ofname)
    print("# N.ch  =", Nch)
    print("# mass1 =", mass1)
    print("# mass2 =", mass2)
    print("# E min =", Emin)
    print("# E del =", Edel)
    print("# E max =", Emax)
    print("# R max =", Rmax)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname, ifengy, ofname, mass1, mass2, Emin, Edel, Emax, Rmax, Nch
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname = ARGV[1].strip()

    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass1'):
                if (Nch != 1):
                    print("\nWARNING: '--mass1' option can be used ONLY for #.ch = 1, pass.\n")
                else:
                    mass1 = float(ARGV[i+1])
            elif (ARGV[i] == '--mass2'):
                if (Nch != 1):
                    print("\nWARNING: '--mass2' option can be used ONLY for #.ch = 1, pass.\n")
                else:
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
            elif (ARGV[i] == '--Nch'):
                Nch = int(ARGV[i+1])
                if (Nch != 1):
                    mass1 = 'From the file'
                    mass2 = 'From the file'
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV[0])
    
    check_args()

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)

    if (argc == 1):
        usage(argv[0])

    set_args(argc, argv)

    if (main() != 0):
        exit("ERROR EXIT.")
