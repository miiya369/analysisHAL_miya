#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Wed Feb  8 16:21:19 JST 2017
### Brief : Calculate the NLO potentials for HAL-potential.

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np

from MiscFuncs.Misc                import *
from MiscFuncs.DataStatistics      import make_mean_err
from MiscFuncs.IO_BinData          import input_bin_data, output_bin_data
from PotentialNLO.CalcPotentialNLO import calc_potential_NLO

Fname_list     = None
Fname_Pot      = None
Fname_Rcorr    = None
Fname_LapRcorr = None

oFname__LO = "/dev/null"
oFname_NLO = "/dev/null"

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [Input Fname list] {--option}\n" % ARGV[0])
    print("option:")
    print("      --fpot   [Potential    file name             ] ... Default ="),; print Fname_Pot
    print("      --fRcorr [R-correlator file name             ] ... Default ="),; print Fname_Rcorr
    print("      --fLapR  [Laplacian R  file name             ] ... Default ="),; print Fname_LapRcorr
    print("      --of_LO  [Output bin-data file name for V__LO]     Default ="),; print oFname__LO
    print("      --of_NLO [Output bin-data file name for V_NLO]     Default ="),; print oFname_NLO
    print("\nNote: [Input Fname list] has not implemented yet, please put dummy.")
    print; quit()

def set_args(ARGC, ARGV):
    global Fname_list, Fname_Pot, Fname_Rcorr, Fname_LapRcorr, oFname__LO, oFname_NLO
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    Fname_list = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--fpot'):
                Fname_Pot = []
                for iarg in range(i+1, ARGC):
                    if (ARGV[iarg][0] != '-'):
                        Fname_Pot.append(ARGV[iarg].strip())
                    else:
                        break
            elif (ARGV[i] == '--fRcorr'):
                Fname_Rcorr = []
                for iarg in range(i+1, ARGC):
                    if (ARGV[iarg][0] != '-'):
                        Fname_Rcorr.append(ARGV[iarg].strip())
                    else:
                        break
            elif (ARGV[i] == '--fLapR'):
                Fname_LapRcorr = []
                for iarg in range(i+1, ARGC):
                    if (ARGV[iarg][0] != '-'):
                        Fname_LapRcorr.append(ARGV[iarg].strip())
                    else:
                        break
            elif (ARGV[i] == '--of_LO'):
                oFname__LO = ARGV[i+1].strip()
            elif (ARGV[i] == '--of_NLO'):
                oFname_NLO = ARGV[i+1].strip()
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile list  ="),; print Fname_list
    print("# ifile Pot   ="),; print Fname_Pot
    print("# ifile Rcor  ="),; print Fname_Rcorr
    print("# ifile LapR  ="),; print Fname_LapRcorr
    print("# ofile V__LO ="),; print oFname__LO
    print("# ofile V_NLO ="),; print oFname_NLO
    print("# =======================")

###### Main part
if __name__ == "__main__":    
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv)
    
    set_args(argc, argv)
    
    tmpNeq1 = len(Fname_Pot)
    tmpNeq2 = len(Fname_Rcorr)
    tmpNeq3 = len(Fname_LapRcorr)
    
    if (not (tmpNeq1 == tmpNeq2 == tmpNeq3)):
        print("\nERROR: #.eq is different, exit.\n"); quit()
    
    check_args()
    
    Neq = tmpNeq1
    
    tmpPot, xData, dummy = input_bin_data(Fname_Pot[0])
    if (tmpPot is xData is dummy is None):
        quit()
    Ndata = len(tmpPot[:, 0])
    Nconf = len(tmpPot[0, :])
    
    Pot      = np.empty((1, 1, Ndata, Neq, Nconf))
    Rcorr    = np.empty((1, 1, Ndata, Neq, Nconf))
    LapRcorr = np.empty((1, 1, Ndata, Neq, Nconf))
    
    for ieq in range(Neq):
        tmpPot, xDataTmp, dummy = input_bin_data(Fname_Pot[ieq])
        for idata in range(Ndata):
            if (xData[idata] != xDataTmp[idata]):
                print("\nERROR: x-data is different, exit.\n"); quit()
            for iconf in range(Nconf):
                Pot[0,0, idata, ieq, iconf] = tmpPot[idata, iconf]
        
        tmpRcorr, xDataTmp, dummy = input_bin_data(Fname_Rcorr[ieq])
        for idata in range(Ndata):
            if (xData[idata] != xDataTmp[idata]):
                print("\nERROR: x-data is different, exit.\n"); quit()
            for iconf in range(Nconf):
                Rcorr[0,0, idata, ieq, iconf] = tmpRcorr[idata, iconf]
        
        tmpLapRcorr, xDataTmp, dummy = input_bin_data(Fname_LapRcorr[ieq])
        for idata in range(Ndata):
            if (xData[idata] != xDataTmp[idata]):
                print("\nERROR: x-data is different, exit.\n"); quit()
            for iconf in range(Nconf):
                LapRcorr[0,0, idata, ieq, iconf] = tmpLapRcorr[idata, iconf]
    
    Pot__LO = np.empty((Ndata, Nconf))
    Pot_NLO = np.empty((Ndata, Nconf))
    
    for iconf in range(Nconf):
        tmpPot__LO, tmpPot_NLO = calc_potential_NLO(Pot     [:,:,:,:,iconf], 
                                                    Rcorr   [:,:,:,:,iconf], 
                                                    LapRcorr[:,:,:,:,iconf])
        if (tmpPot__LO is tmpPot_NLO is None):
            quit()
        for idata in range(Ndata):
            Pot__LO[idata, iconf] = tmpPot__LO[0, 0, idata]
            Pot_NLO[idata, iconf] = tmpPot_NLO[0, 0, idata]
    
    output_bin_data(oFname__LO, Pot__LO, xData)
    output_bin_data(oFname_NLO, Pot_NLO, xData)
    
    for idata in range(Ndata):
        mean0, err0 = make_mean_err(Pot__LO[idata, :])
        mean1, err1 = make_mean_err(Pot_NLO[idata, :])
        print("%lf %e %e %e %e" % (xData[idata], mean0, err0, mean1, err1))
