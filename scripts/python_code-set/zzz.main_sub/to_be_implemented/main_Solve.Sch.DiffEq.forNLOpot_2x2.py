#!/usr/bin/python
# -*- coding: utf-8 -*-

# For LN-SN coupled channel

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np

from MiscFuncs.Misc              import *
from MiscFuncs.DataStatistics    import make_mean_err
from Fitting.IO_Params           import input_params
from Fitting.FitFunctionType     import set_fitfunc_from_fname
from SchDiffeq.SolveDiffEq       import solve_sch_diff_NLO
from Tmatrix.IO_Tmatrix          import output_Tmatrix

hbar_c = 197.327053

iFbase = None
iFengy = None
oFname = '/dev/null'
time   = 13
spin   = 0
mass_N = 1309
mass_L = 1374
mass_S = 1398
Emin   = 0
Edel   = 1
Emax   = 200
Rmax   = 10

lat_a  = 0.0907

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [idir] {--option}\n" % ARGV[0])
    print("option:")
    print("      --ofile [Output file name (For T-matrix)   ] Default ="),; print oFname
    print("      --time  [Time slice to calculate           ] Default ="),; print time
    print("      --spin  [Spin       to calculate           ] Default ="),; print spin
    print("      --Efile [Input file name (For Energy cood.)] Default ="),; print iFengy
    print("      --Emin  [Minimum energy  for output        ] Default ="),; print Emin
    print("      --Edel  [Energy division for output        ] Default ="),; print Edel
    print("      --Emax  [Maximum energy  for output        ] Default ="),; print Emax
    print("      --Rmax  [Maximum  range for output         ] Default ="),; print Rmax
    print("      --lat_a [Lattice spacing (=0 for phys.)    ] Default ="),; print lat_a
    print("\nNote1: The mass, energy and range are should be physical unit (MeV, fm).")
    print; quit()

def set_args(ARGC, ARGV):
    global iFbase, iFengy, oFname, time, spin, Emin, Edel, Emax, Rmax, lat_a
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iFbase = ARGV[1].strip()

    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--ofile'):
                oFname = ARGV[i+1].strip()
            elif (ARGV[i] == '--time'):
                time = int(ARGV[i+1])
            elif (ARGV[i] == '--spin'):
                spin = int(ARGV[i+1])
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
    print("# ifdir  ="),; print iFbase
    print("# if engy="),; print iFengy
    print("# ofile  ="),; print oFname
    print("# time   ="),; print time
    print("# spin   ="),; print spin
    print("# mass N ="),; print mass_N
    print("# mass L ="),; print mass_L
    print("# mass S ="),; print mass_S
    print("# E min  ="),; print Emin
    print("# E del  ="),; print Edel
    print("# E max  ="),; print Emax
    print("# R max  ="),; print Rmax
    print("# lat a  ="),; print lat_a
    print("# =======================")

###### Main part
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv)
    check_args()
    
    # ifiles[LO/NLO, channel(0-3)]
    
    LO__fitfunc = "3G"
    NLO_fitfunc = "2SG"
    
    ifiles = np.array(
        [[iFbase+"/LN-LN_pot_spin%d___LO_t%d_fit_%s"%(spin, time, LO__fitfunc),
          iFbase+"/SN-LN_pot_spin%d___LO_t%d_fit_%s"%(spin, time, LO__fitfunc),
          iFbase+"/SN-LN_pot_spin%d___LO_t%d_fit_%s"%(spin, time, LO__fitfunc),
          iFbase+"/SN-SN_pot_spin%d___LO_t%d_fit_%s"%(spin, time, LO__fitfunc)],
         [iFbase+"/LN-LN_pot_spin%d__NLO_t%d_fit_%s"%(spin, time, NLO_fitfunc),
          iFbase+"/SN-LN_pot_spin%d__NLO_t%d_fit_%s"%(spin, time, NLO_fitfunc),
          iFbase+"/SN-LN_pot_spin%d__NLO_t%d_fit_%s"%(spin, time, NLO_fitfunc),
          iFbase+"/SN-SN_pot_spin%d__NLO_t%d_fit_%s"%(spin, time, NLO_fitfunc)]]
        )
    
    #print ifiles; quit() # For Debug
    
    tmpFuncName, tmpParams = input_params(ifiles[0, 0])
    if (tmpFuncName is tmpParams is None):
        quit()
    
    Nparam = len(tmpParams[:, 0])
    Nconf  = len(tmpParams[0, :])
    Nch    = len(ifiles[0, :]) / 2
    
    # Note: N.param of LO and NLO are assumed to be the same here, 
    #     : but they can be different in this code.
    
    print("#")
    print("# N.conf   ="),; print Nconf
    print("# N.param  ="),; print Nparam
    print("# N.ch     ="),; print Nch
    print("#")
    
### Input the potential func_type and parameters ###
    FitFunc__LO = np.array([[None]*2]*2)
    FitFunc_NLO = np.array([[None]*2]*2)
    Params__LO  = np.empty((Nch, Nch, Nparam, Nconf))
    Params_NLO  = np.empty((Nch, Nch, Nparam, Nconf))
    
    for ich in range(Nch):
        for jch in range(Nch):
            tmpFuncName__LO, tmpParams__LO = input_params(ifiles[0, jch + Nch * ich])
            if (tmpFuncName__LO is tmpParams__LO is None):
                quit()
            tmpFuncName_NLO, tmpParams_NLO = input_params(ifiles[1, jch + Nch * ich])
            if (tmpFuncName_NLO is tmpParams_NLO is None):
                quit()
            
            if (Nparam != len(tmpParams__LO[:, 0]) or Nparam != len(tmpParams_NLO[:, 0])):
                print("\nERROR: #.param is differ, exit.\n"); quit()
            if (Nconf  != len(tmpParams__LO[0, :]) or Nconf  != len(tmpParams_NLO[0, :])):
                print("\nERROR: #.conf is differ, exit.\n"); quit()
            
            FitFunc__LO[ich, jch] = set_fitfunc_from_fname(tmpFuncName__LO)
            FitFunc_NLO[ich, jch] = set_fitfunc_from_fname(tmpFuncName_NLO)
            
            for iparam in range(Nparam):
                for iconf in range(Nconf):
                    Params__LO[ich, jch, iparam, iconf] = tmpParams__LO[iparam, iconf]
                    Params_NLO[ich, jch, iparam, iconf] = tmpParams_NLO[iparam, iconf]
    
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
    Tmat = np.empty((2, 2, N_E, Nconf), dtype=complex)
    mass  = np.array([[mass_N, mass_L], [mass_N, mass_S]])
    
    print("#")
    for iconf in range(Nconf):
        TmpTmat = solve_sch_diff_NLO(FitFunc__LO, Params__LO[:, :, :, iconf], 
                                     FitFunc_NLO, Params_NLO[:, :, :, iconf], 
                                     mass, Edata, Rmax, lat_a)
        for iE in range(N_E):
            for ich in range(2):
                for jch in range(2):
                    Tmat[ich, jch, iE, iconf] = TmpTmat[ich, jch, iE]
        
        print("# Conf = %03d End." % (iconf+1))
    print("#")

### Convert T-matrix ###
    output_Tmatrix(oFname, Tmat, Edata)
