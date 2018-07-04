# -*- coding: utf-8 -*-

"""The module to convert fit parameters for Gaussian-type folding potential."""

def convert_params_gauss(a_Params, a_A, a_RA):
    """
    The function to convert the fit parameters for the Gaussian-type density.
    
    For arguments,
    - a_Params[#.conf, #.param] (2-dim array)
    
    return: r_Params[#.conf, #.param] (2-dim array)
    
    Note: #.param and #.conf are got from a_Params.
    """
    
    from math  import sqrt
    from numpy import empty
    
    l_Nconf  = len(a_Params[:, 0])
    l_Nparam = len(a_Params[0, :])
    
    r_Params = empty((l_Nconf, l_Nparam))
    
    for iconf in range(l_Nconf):
        for iparam in range(0, l_Nparam, 2):
            pn = a_Params[iconf, iparam + 0]
            qn = a_Params[iconf, iparam + 1]
            
            Pn = a_A * pn * pow(qn, 3) / pow(pow(a_RA, 2) + pow(qn, 2), 1.5)
            Qn = sqrt(pow(a_RA, 2) + pow(qn, 2))
            
            r_Params[iconf, iparam + 0] = Pn
            r_Params[iconf, iparam + 1] = Qn
    
    return r_Params
