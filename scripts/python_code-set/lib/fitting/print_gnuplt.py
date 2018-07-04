# -*- coding: utf-8 -*-

"""The module to print fit function form for gnuplot."""

from re import match

def print_gnuplt(a_fname, a_params):
    """
    The function to print the fit function form for gnuplot.
    
    For arguments,
    - a_fname = Simple abridged notation of Fit-function
    - a_params[#.param] (1-dim ndarray)
    
    Note: #.param is got form a_params.
    """
    
    GnuForm = ""
    l_Nparam  = len(a_params)
    
    if   (match('^[1-9]P$', a_fname) is not None): # Power Series
        for ip in range(l_Nparam):
            GnuForm += ("+(%e)*(x**%d)" % (a_params[ip], ip))
    
    elif (match('^Coulomb$', a_fname) is not None): # Coulomb potential
        GnuForm += ("(%e)*(%e)*(1.43996507808)/x)" % (a_params[0], a_params[1]))
    
    elif (match('^[1-9]G$', a_fname) is not None): # Multi-ranges Gaussian
        for ip in range(0, l_Nparam, 2):
            GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[ip], a_params[ip+1]))
    
    elif (match('^1[1-9]G$', a_fname) is not None): # Multi-ranges Gaussian
        for ip in range(0, l_Nparam, 2):
            GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[ip], a_params[ip+1]))
    
    elif (match('^[1-9]SG$', a_fname) is not None): # Multi-ranges shifted-Gaussian
        for ip in range(0, l_Nparam, 3):
            GnuForm += ("+(%e)*exp(-((x-(%e))/(%e))**2)" % 
                        (a_params[ip], a_params[ip+1], a_params[ip+2]))
    
    elif (match('^[1-9]E$', a_fname) is not None): # Multi-ranges Exponential
        for ip in range(0, l_Nparam, 2):
            GnuForm += ("+(%e)*exp(-(%e)*x)" % (a_params[ip], a_params[ip+1]))
    
    elif (match('^[1-9]Y$', a_fname) is not None): # Multi-ranges Yukawa
        for ip in range(0, l_Nparam, 3):
            GnuForm += ("+(%e)*(1-exp(-(%e)*(x**2)))*exp(-(%e)*x)/x" % 
                        (a_params[ip], a_params[ip+1], a_params[ip+2]))
    
    elif (match('^[1-9]Ysq$', a_fname) is not None): # Multi-ranges Yukawa square
        for ip in range(0, l_Nparam, 3):
            GnuForm += ("+(%e)*((1-exp(-(%e)*(x**2)))**2)*exp(-2*(%e)*x)/(x**2)" % 
                        (a_params[ip], a_params[ip+1], a_params[ip+2]))

    elif (match('^[1-9]Ytns$', a_fname) is not None): # Multi-ranges Yukawa tensor
        for ip in range(0, l_Nparam, 3):
            GnuForm += ("+(%e)*((1-exp(-(%e)*(x**2)))**2)*(1+3/((%e)*x)+3/((%e)*x)**2)*exp(-(%e)*x)/x" % 
                        (a_params[ip], a_params[ip+1], a_params[ip+2], a_params[ip+2], a_params[ip+2]))
    
    elif (match('^[1-9]CH$', a_fname) is not None): # Multi-ranges cosh
        for ip in range(0, l_Nparam, 2):
            GnuForm += ("+(%e)*cosh(-(%e)*x)" % (a_params[ip], a_params[ip+1]))
    
    elif (a_fname == "2G1Y"):
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[0], a_params[1]))
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[2], a_params[3]))
        GnuForm += ("+(%e)*(1-exp(-(%e)*(x**2)))*exp(-(%e)*x)/x" % 
                    (a_params[4], a_params[5], a_params[6]))
    
    elif (a_fname == "2G2Y"):
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[0], a_params[1]))
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[2], a_params[3]))
        GnuForm += ("+(%e)*(1-exp(-(%e)*(x**2)))*exp(-(%e)*x)/x" % (a_params[4], a_params[5], a_params[6]))
        GnuForm += ("+(%e)*(1-exp(-(%e)*(x**2)))*exp(-(%e)*x)/x" % (a_params[7], a_params[8], a_params[9]))
    
    elif (a_fname == "2G1Ysq"):
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[0], a_params[1]))
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[2], a_params[3]))
        GnuForm += ("+(%e)*((1-exp(-(%e)*(x**2)))**2)*exp(-2*(%e)*x)/(x**2)" % 
                    (a_params[4], a_params[5], a_params[6]))

    elif (a_fname == "3G1Ysq"):
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[0], a_params[1]))
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[2], a_params[3]))
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[4], a_params[5]))
        GnuForm += ("+(%e)*((1-exp(-(%e)*(x**2)))**2)*exp(-2*(%e)*x)/(x**2)" % 
                    (a_params[6], a_params[7], a_params[8]))
    
    elif (a_fname == "2G1Y1Ysq"):
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[0], a_params[1]))
        GnuForm += ("+(%e)*exp(-(x/(%e))**2)" % (a_params[2], a_params[3]))
        GnuForm += ("+(%e)*(1-exp(-(%e)*(x**2)))*exp(-(%e)*x)/x" % (a_params[4], a_params[5], a_params[6]))
        GnuForm += ("+(%e)*((1-exp(-(%e)*(x**2)))**2)*exp(-2*(%e)*x)/(x**2)" % 
                    (a_params[7], a_params[8], a_params[9]))
    
    else:
        print("# Warning: Invalid (or Have not implemented) fit function type: %s" % a_fname)
        GnuForm += "None"
    
    print("# For Gnuplot: %s" % GnuForm)
