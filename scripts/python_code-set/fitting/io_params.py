# -*- coding: utf-8 -*-

"""The module to input/output fit parameters."""

from numpy import empty, append, inf, loadtxt

def input_init_params(a_iFname):
    """
    The function to input the initial-fit parameters.
    
    FORMAT: [parameter] [lower bound (if blank, -inf)] [upper bound (if blank, inf)]
    FORMAT: #.row corresponds to #.parameter
    FORMAT: A blank row and the row starting from '#' will be skippled
    
    return: (Initial Params[#.param], 
             Lower bounds  [#.param], 
             Upper bounds  [#.param])
    """
    
    r_in_Params = empty(0)
    r_Upper_b   = empty(0)
    r_Lower_b   = empty(0)
    
    for Line in open(a_iFname, 'r'):
        line = Line.split()
        if (len(line) == 0 or line[0][0].strip() == '#'):
            continue
        
        if (len(line) == 1 or line[1][0].strip() == '#'):
            r_in_Params = append(r_in_Params, float(line[0].strip()))
            r_Lower_b   = append(r_Lower_b  , -inf)
            r_Upper_b   = append(r_Upper_b  ,  inf)
        elif (len(line) == 2 or line[2][0].strip() == '#'):
            r_in_Params = append(r_in_Params, float(line[0].strip()))
            r_Lower_b   = append(r_Lower_b  , float(line[1].strip()))
            r_Upper_b   = append(r_Upper_b  , inf)
        elif (len(line) == 3 or line[3][0].strip() == '#'):
            r_in_Params = append(r_in_Params, float(line[0].strip()))
            r_Lower_b   = append(r_Lower_b  , float(line[1].strip()))
            r_Upper_b   = append(r_Upper_b  , float(line[2].strip()))
        else:
            print("\nERROR: Too many rows in an initial fit-parameter file, exit.\n")
            return (None, None, None)
    
    return (r_in_Params, r_Lower_b, r_Upper_b)

def output_params(a_oFname, a_FuncName, a_Params, verbose_flg = True):
    """
    The function to output the fit parameters.
    
    FORMAT: First line = '#FitFunction: [Fit function name]'
    FORMAT: From secound line, #.row = #.conf, #.column = #.param
    
    For arguments,
    - a_FuncNname = Simple abridged notation of Fit-function
    - a_Params[#.conf, #.param] (2-dim ndarray)
    
    Note: #.conf and #.param are got from a_Params.
    """
    
    l_Nconf  = len(a_Params[:,0])
    l_Nparam = len(a_Params[0,:])
    
    with open(a_oFname, 'w') as ofile:
        ofile.write("#FitFunction: %s\n" % a_FuncName)
        
        for iconf in range(l_Nconf):
            ofile.write("%1.16e "*len(a_Params[iconf,:]) % (*a_Params[iconf,:],) + "\n")
    
    if (verbose_flg):
        print("# Successful to output fit parameters.")

def input_params(a_iFname, verbose_flg = True):
    """
    The function to input the fit parameters.
    
    FORMAT: First line = '#FitFunction: [Fit function name]'
    FORMAT: From secound line, #.row = #.conf, #.column = #.param
    
    return: (FitFunc Name, Params[#.conf, #.param])
    """
    
    from fitting.fitfunc_type import get_Nparam_from_fname
    
    for Line in open(a_iFname, 'r'):
        if (Line.split()[0].strip() != '#FitFunction:'):
            print("\nERROR: Invalid fit-parameter file, exit.\n")
            return (None, None)
        r_FuncName = Line.split()[1].strip()
        break
    
    if (get_Nparam_from_fname(r_FuncName) is None):
        return (None, None)
    
    l_Nparam = get_Nparam_from_fname(r_FuncName)
    r_Params = loadtxt(a_iFname)
    if (len(r_Params[0,:]) != l_Nparam):
        print("\nERROR: Unexpected #.row or #.column in the parameter file, exit.\n")
        return (None, None)
    
    if (verbose_flg):
        print("# Successful to input fit parameters.")
        print("# N.conf   = %d" % len(r_Params[:,0]))
        print("# N.param  = %d" % l_Nparam)
        print("# FuncName = %s" % r_FuncName)
    
    return (r_FuncName, r_Params)
