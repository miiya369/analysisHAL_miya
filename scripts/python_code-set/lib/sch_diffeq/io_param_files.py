# -*- coding: utf-8 -*-

"""The module to input NxN fit parameters lists."""

from numpy             import array
from fitting.io_params import input_params

def input_NxN_params(a_iFname, verbose_flg = True):
    """
    The function to input the initial-fit parameters.
    
    FORMAT: First  line = '#.ch: [Number of channel]'
    FORMAT: Second line = [mass of ch.1-1] [mass of ch.1-2] ... [mass of ch.N-1] [mass of ch.N-2]
    FORMAT: From third line, [a name of parameter file] (for each line)
    FORMAT: The order of file name corresponds to jch+Nch*ich for params[ich, jch]
    FORMAT: A blank row and the row starting from '#' will be skippled
    
    return: (FuncName[#.ch, #.ch],
             Params  [#.ch, #.ch, #.conf, #.param],
             mass    [#.ch, 2])
    """    
    with open(a_iFname, 'r') as ifile:
        lines = ifile.readlines()
    
    if (lines[0].split()[0].strip() != '#.ch:'):
        print("\nERROR: Invalid NxN fit-parameters file (error at first line), exit.\n")
        return (None, None, None)
    
    l_Nch = int(lines[0].split()[1].strip())
    
    if (len(lines[1].split()) != l_Nch*2):
        print("\nERROR: Invalid NxN fit-parameters file (error at second line), exit.\n")
        return (None, None, None)
    
    r_mass = array([[float(lines[1].split()[j+2*i].strip()) 
                     for j in range(2)] for i in range(l_Nch)])
    
    l_fparams = []; tmp_count = 0
    for i in range(2, len(lines)):
        line = lines[i].split()
        if (len(line) == 0 or line[0][0].strip() == '#'):
            continue
        l_fparams.append(line[0])
        tmp_count += 1
        if (tmp_count == l_Nch**2):
            break
    
    if (len(l_fparams) != l_Nch**2):
        print("\nERROR: Invalid NxN fit-parameters file (error at below third line), exit.\n")
        return (None, None, None)
    
    r_FuncName = array([[input_params(l_fparams[jch+l_Nch*ich], verbose_flg=False)[0] 
                         for jch in range(l_Nch)] for ich in range(l_Nch)])
    
    r_Params   = array([[input_params(l_fparams[jch+l_Nch*ich], verbose_flg=False)[1] 
                         for jch in range(l_Nch)] for ich in range(l_Nch)])
    
    tmp_Nconf = len(r_Params[0,0,:,0])
    for ich in range(l_Nch):
        for jch in range(l_Nch):
            if (tmp_Nconf != len(r_Params[ich,jch,:,0])):
                print("\nERROR: Invalid NxN fit-parameters file (different #.conf), exit.\n")
                return (None, None, None)
    
    if (verbose_flg):
        print("# Successful to input NxN fit parameters.")
        print("# N.ch     =", l_Nch)
        print("# N.conf   =", tmp_Nconf)
        for ich in range(l_Nch):
            print("# mass[%2d] =" % ich, r_mass[ich,:])
        print("# Fit function (%dx%d):" % (l_Nch, l_Nch))
        for ich in range(l_Nch):
            print("#", "%9s"*l_Nch % (*r_FuncName[ich,:],))
    
    return (r_FuncName, r_Params, r_mass)
