# -*- coding: utf-8 -*-

"""The module to convert fit parameters for Gaussian-type folding potential."""

def convert_params_gauss(a_Params, a_A, a_RA):
    """
    The function to convert the fit parameters for the Gaussian-type density.
    
    For arguments,
    - a_Params[#.param, #.conf] (2-dim array)
    
    return: r_Params[#.param, #.conf] (2-dim array)
    
    Note: #.param and #.data are got from a_Params.
    """
    
    from math  import sqrt
    from numpy import empty
    
    l_Nparam = len(a_Params[:, 0])
    l_Nconf  = len(a_Params[0, :])
    
    r_Params = empty((l_Nparam, l_Nconf))
    
    for iconf in range(l_Nconf):
        for iparam in range(0, l_Nparam, 2):
            pn = a_Params[iparam + 0, iconf]
            qn = a_Params[iparam + 1, iconf]
            
            Pn = a_A * pn * pow(qn, 3) / pow(pow(a_RA, 2) + pow(qn, 2), 1.5)
            Qn = sqrt(pow(a_RA, 2) + pow(qn, 2))
            
            r_Params[iparam + 0, iconf] = Pn
            r_Params[iparam + 1, iconf] = Qn
    
    return r_Params
