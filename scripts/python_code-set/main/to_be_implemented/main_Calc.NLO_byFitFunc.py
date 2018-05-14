#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from MiscFuncs.Misc           import *
from MiscFuncs.DataStatistics import make_mean_err
from Fitting.IO_Params        import input_params
from Fitting.FitFunctionType  import set_fitfunc_from_fname
from MiscFuncs.IO_BinData     import input_bin_data, output_bin_data

iFname_fit  = None
iFname_bin1 = None
iFname_bin2 = None
iFname_bin3 = None
iFname_bin4 = None

oFname__LO = "/dev/null"
oFname_NLO = "/dev/null"

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [fit params (V_wall-V_exp)] [bin data 1] [bin data 2] [bin data 3] [bin data 4] {--option}\n" % ARGV[0])
    print("   - bin-data 1 = V_wall")
    print("   - bin-data 2 = V_exp")
    print("   - bin-data 3 = Lap.R/R (wall)")
    print("   - bin-data 4 = Lap.R/R (exp)\n")
    print("option:")
    print("      --of_LO  [Output bin-data file name for V__LO] Default ="),; print oFname__LO
    print("      --of_NLO [Output bin-data file name for V_NLO] Default ="),; print oFname_NLO
    print; quit()

def set_args(ARGC, ARGV):
    global iFname_fit, iFname_bin1, iFname_bin2, iFname_bin3, iFname_bin4, oFname__LO, oFname_NLO
    
    if (ARGV[1][0] == '-' or ARGV[2][0] == '-' or ARGV[3][0] == '-' or ARGV[4][0] == '-' or ARGV[5][0] == '-'):
        usage(ARGV)
    
    iFname_fit  = ARGV[1].strip()
    iFname_bin1 = ARGV[2].strip()
    iFname_bin2 = ARGV[3].strip()
    iFname_bin3 = ARGV[4].strip()
    iFname_bin4 = ARGV[5].strip()
    
    for i in range(6, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--of_LO'):
                oFname__LO = ARGV[i+1].strip()
            elif (ARGV[i] == '--of_NLO'):
                oFname_NLO = ARGV[i+1].strip()
            else:
                print("\nERROR: Invalid option '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile fit   ="),; print iFname_fit
    print("# ifile bin 1 ="),; print iFname_bin1
    print("# ifile bin 2 ="),; print iFname_bin2
    print("# ifile bin 3 ="),; print iFname_bin3
    print("# ifile bin 4 ="),; print iFname_bin4
    print("# ofile V__LO ="),; print oFname__LO
    print("# ofile V_NLO ="),; print oFname_NLO
    print("# =======================")

###### Main part
if __name__ == "__main__":
    import numpy as np
    
    argv = sys.argv; argc = len(argv)
    
    if (argc < 6):
        usage(argv)
    
    set_args(argc, argv)
    check_args()
    
    Fname, Params = input_params(iFname_fit)
    if (Fname is Params is None):
        quit()
    
    Nconf   = len(Params[0, :])
    Nparam  = len(Params[:, 0])
    FitFunc = set_fitfunc_from_fname(Fname)
    
    yData1, xData1, Dummy = input_bin_data(iFname_bin1)
    if (xData1 is yData1 is Dummy is None):
        quit()
    if (Nconf != len(yData1[0, :])):
        print("\nERROR: #.conf is differ, exit.\n"); quit()
    
    yData2, xData2, Dummy = input_bin_data(iFname_bin2)
    if (xData2 is yData2 is Dummy is None):
        quit()
    if (Nconf != len(yData2[0, :])):
        print("\nERROR: #.conf is differ, exit.\n"); quit()
    
    yData3, xData3, Dummy = input_bin_data(iFname_bin3)
    if (xData3 is yData3 is Dummy is None):
        quit()
    if (Nconf != len(yData3[0, :])):
        print("\nERROR: #.conf is differ, exit.\n"); quit()
    
    yData4, xData4, Dummy = input_bin_data(iFname_bin4)
    if (xData4 is yData4 is Dummy is None):
        quit()
    if (Nconf != len(yData4[0, :])):
        print("\nERROR: #.conf is differ, exit.\n"); quit()
    
    if (len(xData1) != len(xData2) or 
        len(xData1) != len(xData3) or 
        len(xData1) != len(xData4) or 
        len(xData2) != len(xData3) or 
        len(xData2) != len(xData4) or 
        len(xData3) != len(xData4)):
        print("\nERROR: #.data is differ, exit.\n"); quit()
    
    Ndata = len(xData1)
    
    V__LO = np.empty((Ndata, Nconf))
    V_NLO = np.empty((Ndata, Nconf))
    
    tmpV_LO1 = np.empty(Nconf)
    tmpV_LO2 = np.empty(Nconf)
    
    for i in range(Ndata):
        if (xData1[i] != xData2[i] or 
            xData1[i] != xData3[i] or 
            xData1[i] != xData4[i] or 
            xData2[i] != xData3[i] or 
            xData2[i] != xData4[i] or 
            xData3[i] != xData4[i]):
            print("\nERROR: x-data is differ, exit.\n"); quit()
        
        for iconf in range(Nconf):
            V_NLO[i, iconf] = FitFunc(xData1[i], *Params[:, iconf]) / (yData3[i, iconf] - yData4[i, iconf])
            
            tmpV_LO1[iconf] = yData1[i, iconf] - yData3[i, iconf] * V_NLO[i, iconf]
            tmpV_LO2[iconf] = yData2[i, iconf] - yData4[i, iconf] * V_NLO[i, iconf]
            
            # Selection:
            V__LO[i, iconf] =  tmpV_LO1[iconf]
            #V__LO[i, iconf] =  tmpV_LO2[iconf]
            #V__LO[i, iconf] = (tmpV_LO1[iconf] + tmpV_LO2[iconf]) / 2.0
        
        m__LO, e__LO = make_mean_err(V__LO[i, :])
        m_NLO, e_NLO = make_mean_err(V_NLO[i, :])
        # For Debug:
        #m__LO, e__LO = make_mean_err(tmpV_LO1)
        #m_NLO, e_NLO = make_mean_err(tmpV_LO2)
        print("%lf %1.16e %1.16e %1.16e %1.16e" % (xData1[i], m__LO, e__LO, m_NLO, e_NLO))
    
    output_bin_data(oFname__LO, V__LO, xData1)
    output_bin_data(oFname_NLO, V_NLO, xData1)
