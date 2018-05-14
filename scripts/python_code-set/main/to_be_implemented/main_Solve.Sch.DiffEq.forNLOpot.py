#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np

from MiscFuncs.Misc              import *
from MiscFuncs.DataStatistics    import make_mean_err
from Fitting.IO_Params           import input_params
from Fitting.FitFunctionType     import set_fitfunc_from_fname
from SchDiffeq.SolveDiffEq       import solve_sch_diff_NLO
from Tmatrix.ConvertMat          import convert_TtoS
from Tmatrix.CalcPhaseShift_Smat import calc_phase_Nch1

hbar_c = 197.327053

iF__LO = None
iF_NLO = None
iFengy = None
oFname = '/dev/null'
mass1  = 1308.924032
mass2  = 1374.487402
Emin   = 0
Edel   = 1
Emax   = 200
Rmax   = 10

lat_a  = 0.0907

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [fparam (V_LO [MeV])] [fparam (V_NLO [MeV])] {--option}\n" % ARGV[0])
    print("option:")
    print("      --ofile [Output file name (For T-matrix)   ] Default ="),; print oFname
    print("      --mass1 [Mass of 1st baryon                ] Default ="),; print mass1
    print("      --mass2 [Mass of 2nd baryon                ] Default ="),; print mass2
    print("      --Efile [Input file name (For Energy cood.)] Default ="),; print iFengy
    print("      --Emin  [Minimum energy  for output        ] Default ="),; print Emin
    print("      --Edel  [Energy division for output        ] Default ="),; print Edel
    print("      --Emax  [Maximum energy  for output        ] Default ="),; print Emax
    print("      --Rmax  [Maximum  range for output         ] Default ="),; print Rmax
    print("      --lat_a [Lattice spacing                   ] Default ="),; print lat_a
    print("\nNote: The mass, energy and range are should be physical unit (MeV, fm).")
    print; quit()

def set_args(ARGC, ARGV):
    global iF__LO, iF_NLO, iFengy, oFname, mass1, mass2, Emin, Edel, Emax, Rmax, lat_a
    
    if (ARGV[1][0] == '-' or ARGV[2][0] == '-'):
        usage(ARGV)
    
    iF__LO = ARGV[1].strip()
    iF_NLO = ARGV[2].strip()

    for i in range(3, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass1'):
                mass1 = float(ARGV[i+1])
            elif (ARGV[i] == '--mass2'):
                mass2 = float(ARGV[i+1])
            elif (ARGV[i] == '--ofile'):
                oFname = ARGV[i+1].strip()
            elif (ARGV[i] == '--Efile'):
                iFengy = ARGV[i+1].strip()
                Emin   = 'From the file'
                Edel   = 'From the file'
                Emax   = 'From the file'
            elif (ARGV[i] == '--Emin'):
                if (iFengy is not None):
                    print("\nWARNING: '--Emin' option is ignored " +
                          "when '--Efile' option was specified.\n")
                else:
                    Emin = float(ARGV[i+1])
            elif (ARGV[i] == '--Edel'):
                if (iFengy is not None):
                    print("\nWARNING: '--Edel' option is ignored " +
                          "when '--Efile' option was specified.\n")
                else:
                    Edel = float(ARGV[i+1])
            elif (ARGV[i] == '--Emax'):
                if (iFengy is not None):
                    print("\nWARNING: '--Emax' option is ignored " +
                          "when '--Efile' option was specified.\n")
                else:
                    Emax = float(ARGV[i+1])
            elif (ARGV[i] == '--Rmax'):
                Rmax = float(ARGV[i+1])
            elif (ARGV[i] == '--lat_a'):
                lat_a = float(ARGV[i+1])
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# if  LO ="),; print iF__LO
    print("# if NLO ="),; print iF_NLO
    print("# ifengy ="),; print iFengy
    print("# ofile  ="),; print oFname
    print("# mass1  ="),; print mass1
    print("# mass2  ="),; print mass2
    print("# E min  ="),; print Emin
    print("# E del  ="),; print Edel
    print("# E max  ="),; print Emax
    print("# R max  ="),; print Rmax
    print("# lat a  ="),; print lat_a
    print("# =======================")

###### Main part
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc < 3):
        usage(argv)
    
    set_args(argc, argv)
    check_args()
    
### Input the potential func_type and parameters ###
    FuncName__LO, tmpParams__LO = input_params(iF__LO)
    if (FuncName__LO is tmpParams__LO is None):
        quit()
    FuncName_NLO, tmpParams_NLO = input_params(iF_NLO)
    if (FuncName_NLO is tmpParams_NLO is None):
        quit()
    
    if (len(tmpParams__LO[0, :]) != len(tmpParams_NLO[0, :])):
        print("\nERROR: #.conf is differ, exit.\n"); quit()
    
    Nconf = len(tmpParams__LO[0, :])
    mass  = np.array([[mass1, mass2]])
    
    Params__LO = np.empty((1, 1, len(tmpParams__LO[:, 0]), Nconf))
    Params_NLO = np.empty((1, 1, len(tmpParams_NLO[:, 0]), Nconf))
    for iconf in range(Nconf):
        for iparam in range(len(tmpParams__LO[:, 0])):
            Params__LO[0, 0, iparam, iconf] = tmpParams__LO[iparam, iconf]
        for iparam in range(len(tmpParams_NLO[:, 0])):
            Params_NLO[0, 0, iparam, iconf] = tmpParams_NLO[iparam, iconf]
    
### Set the energies to calculate
    if (iFengy is None):
        N_E   = int((Emax - Emin) / Edel)
        Edata = np.empty(N_E)
        for i in range(N_E):
            Edata[i] = Emin + i*Edel
    else:
        ifile = open(iFengy, 'r')
        Lines = ifile.readlines()
        ifile.close()

        N_E   = len(Lines)
        Edata = np.empty(N_E)
        for i in range(N_E):
            Edata[i] = float(Lines[i].split()[0].strip())

    #print Edata; quit() # For Debug
    
### T-matrix calculation ###
    Tmat = np.empty((1, 1, N_E, Nconf), dtype=complex)
    
    print("#")
    for iconf in range(Nconf):        
        TmpTmat = solve_sch_diff_NLO(np.array([[set_fitfunc_from_fname(FuncName__LO)]]), 
                                     Params__LO[:, :, :, iconf], 
                                     np.array([[set_fitfunc_from_fname(FuncName_NLO)]]), 
                                     Params_NLO[:, :, :, iconf], 
                                     mass, Edata, Rmax, lat_a)
        
        for iE in range(N_E):
            Tmat[0, 0, iE, iconf] = TmpTmat[0, 0, iE]
        
        print("# Conf = %03d End." % (iconf+1))
    print("#")

### Output the Phase shift & Calculation of the Scattering length ###
    convert_TtoS(Tmat)

    tmp_m      = (mass1*mass2)/(mass1+mass2)
    PhaseShift = np.empty(Nconf)
    ScatLength = np.empty(Nconf)

    print("# E, phs, phs_e, scatt.len, scatt.len_e")
    for iE in range(N_E):
        if (Edata[iE] == 0.0):
            Edata[iE] = 10e-10
        for iconf in range(Nconf):
            PhaseShift[iconf] = calc_phase_Nch1(Tmat[0, 0, iE, iconf])[0]
            ScatLength[iconf] = (np.tan(PhaseShift[iconf] * np.pi/180) /
                                 np.sqrt(2.0 * tmp_m * Edata[iE]) * hbar_c)
        if (Edata[iE] == 10e-10):
            Edata[iE] = 0.0
        
        mean1, err1 = make_mean_err(PhaseShift)
        mean2, err2 = make_mean_err(ScatLength)
        print("%lf %e %e %e %e" % (Edata[iE], mean1, err1, mean2, err2))

    ScatLength_i = np.empty(Nconf)
    tmp_E        = 1e-10
    for iconf in range(Nconf):        
        TmpTmat = solve_sch_diff_NLO(np.array([[set_fitfunc_from_fname(FuncName__LO)]]), 
                                     Params__LO[:, :, :, iconf], 
                                     np.array([[set_fitfunc_from_fname(FuncName_NLO)]]), 
                                     Params_NLO[:, :, :, iconf], 
                                     mass, np.array([tmp_E]), Rmax, lat_a)[0, 0, 0]
        
        # Scattering length [fm] := T(p)/p (p->0)
        ScatLength  [iconf] = TmpTmat.real / np.sqrt(2.0 * tmp_m * tmp_E) * hbar_c
        ScatLength_i[iconf] = TmpTmat.imag / np.sqrt(2.0 * tmp_m * tmp_E) * hbar_c

    m_r, e_r = make_mean_err(ScatLength  )
    m_i, e_i = make_mean_err(ScatLength_i)
    print("#\n# Scattering length [fm] = %lf(%lf) + %lf(%lf) i" % (m_r, e_r, m_i, e_i))
