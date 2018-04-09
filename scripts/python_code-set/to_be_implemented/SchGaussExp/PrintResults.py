# -*- coding: utf-8 -*-

"""The module to output the results for Gauss expansion method."""

from MiscFuncs.DataStatistics import make_mean_err

def print_wave(a_EigVec, a_Range, a_max_r, a_div_r = 0.01):
    """
    The function to output the wave functions for Gauss expansion method.
    
    For arguments,
    - a_EigVec[#.base, #.return, #.conf] (3-dim array)
    - a_Range [#.base]                   (1-dim array)
    
    Note: #.base, #.return and #.conf are got from a_EigVec.
    """
    
    from numpy          import sum, exp, zeros
    from MiscFuncs.Misc import frange
    
    OUT_WAVE_FLG = True # For Debug
    
    l_Nbase = len(a_EigVec[:, 0, 0])
    l_Nret  = len(a_EigVec[0, :, 0])
    l_Nconf = len(a_EigVec[0, 0, :])
    
    print("# Range     ="),
    for i in range(l_Nbase):
        print("%1.16e" % a_Range[i]),
    print("")
    for iret in range(l_Nret):
        print("#\n# Return %d" % iret)
        for iconf in range(l_Nconf):
            print("# Coeff %3d =" % iconf),
            for i in range(l_Nbase):
                print("%1.16e" % a_EigVec[i, iret, iconf]),
            print("")
    
    if (OUT_WAVE_FLG): # For Debug
        print("#")
        tmp_d = zeros(l_Nconf)
        for r in frange(0.0, a_max_r, a_div_r):
            print("%lf" % r),
            for iret in range(l_Nret):
                for iconf in range(l_Nconf):
                    tmp_d[iconf] = sum(a_EigVec[:, iret, iconf] * exp(-pow(r/a_Range, 2)))
                mean, err = make_mean_err(tmp_d)
                print("%e %e" % (mean, err)),
            print("")
        del tmp_d

def print_results(a_EigVal, a_EigVec, a_Range, a_max_r, a_div_r = 0.01, a_EigVal_Only = True):
    """
    The function to output the results for Gauss expansion method.
    
    For arguments,
    - a_EigVal[#.return, #.conf]           (2-dim array)
    - a_EigVec[#.base  , #.return, #.conf] (3-dim array)
    - a_Range [#.base]                     (1-dim array)
    
    Note: #.return is got from a_EigVal.
    """
    
    from numpy import zeros, append
    
    l_Nconf = len(a_EigVal[0, :])
    l_Nret  = len(a_EigVal[:, 0])
    
    if (not a_EigVal_Only):
        print("#")
        print_wave(a_EigVec, a_Range, a_max_r, a_div_r)
    
    print("#")
    for iret in range(l_Nret):
        tmpRes = zeros(0)
        for iconf in range(l_Nconf):
            #if (a_EigVal[iret, iconf] < -0.1):
            tmpRes = append(tmpRes, a_EigVal[iret, iconf])
        
        if (len(tmpRes) == 0):
            mean = 0.0
            err  = 0.0
        elif (len(tmpRes) == 1):
            mean = tmpRes[0]
            err  = 0.0
        else:
            mean, err = make_mean_err(tmpRes)#, is_JKdata = False)
        
        print("# Eigen energy[%d] = %lf +/- %lf" % (iret, mean, err))
        
        #print("# Eigen energy[%d] = " % iret), # For All data
        #for iconf in range(l_Nconf):
        #    if(a_EigVal[iret, iconf] < 0):
        #        print(" %lf" % a_EigVal[iret, iconf]),
        #    else:
        #        print(" 0.0"),
        #print
