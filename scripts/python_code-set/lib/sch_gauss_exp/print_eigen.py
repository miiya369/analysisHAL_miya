# -*- coding: utf-8 -*-

"""The module to print the eigen values or eigen vectors for Gauss expansion method."""

from __future__        import print_function
from common.statistics import make_mean_err

def print_wave(a_evec, a_range, a_max_r, a_div_r = 0.01, owave_flg = True):
    """
    The function to output the wave functions for Gauss expansion method.
    
    For arguments,
    - a_evec [#.conf, #.base, #.state] (3-dim ndarray)
    - a_range[#.base]                  (1-dim ndarray)
    
    Note: #.base, #.state and #.conf are got from a_evec.
    """    
    from numpy       import sum, exp, array
    from common.misc import frange
    
    l_Nconf = len(a_evec[:,0,0])
    l_Nbase = len(a_evec[0,:,0])
    l_Nstat = len(a_evec[0,0,:])
    
    print("# range =", ("%1.16e "*len(a_range)) % (*a_range,))
    
    for istat in range(l_Nstat):
        print("#\n# State %03d:" % istat)
        for iconf in range(l_Nconf):
            print("# Coeffs %03d =" % iconf, ("%1.16e "*l_Nbase) % (*a_evec[iconf,:,istat],))
    
    if (owave_flg):
        print("#")
        for r in frange(0.0, a_max_r, a_div_r):
            tmp1 = array([[sum(a_evec[iconf,:,istat] * exp(-pow(r/a_range, 2))) 
                           for iconf in range(l_Nconf)] for istat in range(l_Nstat)])
            tmp2 = array([make_mean_err(tmp1[istat,:])  for istat in range(l_Nstat)])
            print("%lf" % r, ("%e %e "*len(tmp2[:,0])) % (*tmp2.flatten(),))

def print_eigen(a_eval, a_evec, a_range, a_max_r, a_div_r = 0.01, a_eval_only = True):
    """
    The function to output the results for Gauss expansion method.
    
    For arguments,
    - a_eval [#.conf, #.state]         (2-dim ndarray)
    - a_evec [#.conf, #.base, #.state] (3-dim ndarray)
    - a_range[#.base]                  (1-dim ndarray)
    
    Note: #.state is got from a_eval.
    """    
    if (not a_eval_only):
        print("#")
        print_wave(a_evec, a_range, a_max_r, a_div_r)
    
    l_Nconf = len(a_eval[:,0])
    l_Nstat = len(a_eval[0,:])
    
    print("#")
    for istat in range(l_Nstat):
        print("# Eigen energy[%02d] = %lf +/- %lf" % (istat, *make_mean_err(a_eval[:,istat])))
        ### For All data
        #print("# Eigen energy[%02d] =" % istat, ("%lf "*l_Nconf) % (*a_eval[:,istat],))
