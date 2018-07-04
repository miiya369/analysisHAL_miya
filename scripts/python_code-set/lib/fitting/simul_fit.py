# -*- coding: utf-8 -*-

"""The module to verify the fit functions type for the simultaniously fitting."""

from numpy import exp

def fitfunc_3G_Simul(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9):
    """3-ranges Gaussian for the simultaniously fitting"""

    return (a_p1 * exp(-(a_x/a_p7)**2) + a_p2 * exp(-(a_x/a_p8)**2) +
            a_p3 * exp(-(a_x/a_p9)**2) + 
            a_p4 * exp(-((a_x-100)/a_p7)**2) + a_p5 * exp(-((a_x-100)/a_p8)**2) + 
            a_p6 * exp(-((a_x-100)/a_p9)**2))

def fitfunc_4G_Simul(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9, a_p10, a_p11, a_p12):
    """4-ranges Gaussian for the simultaniously fitting"""

    return (a_p1 * exp(-(a_x/a_p9 )**2) + a_p2 * exp(-(a_x/a_p10)**2) +
            a_p3 * exp(-(a_x/a_p11)**2) + a_p4 * exp(-(a_x/a_p12)**2) +
            a_p5 * exp(-((a_x-100)/a_p9 )**2) + a_p6 * exp(-((a_x-100)/a_p10)**2) + 
            a_p7 * exp(-((a_x-100)/a_p11)**2) + a_p8 * exp(-((a_x-100)/a_p12)**2))

def fitfunc_2G1Ysq_Simul(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9, a_p10):
    """2-ranges Gaussian + one Yukawa square for the simultaniously fitting (Assume that r != 0)"""
    
    return (a_p1 * exp(-pow(a_x/a_p7, 2)) + a_p2 * exp(-pow(a_x/a_p8, 2)) +
            a_p3 * (1.0 - exp(-a_p9 * a_x**2))**2 * exp(-2.0 * a_p10 * a_x) / a_x**2 +
            a_p4 * exp(-pow((a_x-100)/a_p7, 2)) + a_p5 * exp(-pow((a_x-100)/a_p8, 2)) +
            a_p6 * (1.0 - exp(-a_p9 * (a_x-100)**2))**2 * exp(-2.0 * a_p10 * (a_x-100)) / (a_x-100)**2)

def fitfunc_3G1Ysq_Simul(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9, a_p10, a_p11, a_p12, a_p13):
    """3-ranges Gaussian + one Yukawa square for the simultaniously fitting (Assume that r != 0)"""
    
    return (a_p1 * exp(-pow(a_x/a_p9, 2)) + a_p2 * exp(-pow(a_x/a_p10, 2)) +
            a_p3 * exp(-pow(a_x/a_p11, 2)) +
            a_p4 * (1.0 - exp(-a_p12 * a_x**2))**2 * exp(-2.0 * a_p13 * a_x) / a_x**2 +
            a_p5 * exp(-pow((a_x-100)/a_p9, 2)) + a_p6 * exp(-pow((a_x-100)/a_p10, 2)) +
            a_p7 * exp(-pow((a_x-100)/a_p11, 2)) +
            a_p8 * (1.0 - exp(-a_p12 * (a_x-100)**2))**2 * exp(-2.0 * a_p13 * (a_x-100)) / (a_x-100)**2)

def get_Nparam_from_fname_simul(a_fname):
    """
    The function to get the #.param of various fit functions for the simultaniously fitting.
    (Fit fname ==> #.param)
    
    return: #.param
    
    Note: The definition of the fname is noted in ./Definition_Ftype.txt
    """
    
    if   (a_fname == '3G_Simul'):
        return 9
    elif (a_fname == '2G1Ysq_Simul'):
        return 10
    elif (a_fname == '4G_Simul'):
        return 12
    elif (a_fname == '3G1Ysq_Simul'):
        return 13
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_fname)
        return None

def set_fitfunc_from_fname_simul(a_fname):
    """
    The function to get the fit function for the simultaniously fitting.
    (Fit fname ==> Fit function)
    
    return: Fit function
    
    Note: The definition of the fname is noted in ./Definition_Ftype.txt
    """
    
    if   (a_fname == '3G_Simul'):
        return fitfunc_3G_Simul
    elif (a_fname == '4G_Simul'):
        return fitfunc_4G_Simul
    elif (a_fname == '2G1Ysq_Simul'):
        return fitfunc_2G1Ysq_Simul
    elif (a_fname == '3G1Ysq_Simul'):
        return fitfunc_3G1Ysq_Simul
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_fname)
        return None

def convert_simul_fitparam(a_ps, a_fname, a_idx):
    """
    The function to convert the parameters (simultanious --> divided).
    """
    
    if   (a_fname == '3G_Simul'):
        return ("3G", (a_ps[a_idx*3], a_ps[6], a_ps[a_idx*3+1], a_ps[7], a_ps[a_idx*3+2], a_ps[8]))
    elif (a_fname == '4G_Simul'):
        return ("4G", (a_ps[a_idx*4], a_ps[8], a_ps[a_idx*4+1], a_ps[9], a_ps[a_idx*4+2], a_ps[10], 
                       a_ps[a_idx*4+3], a_ps[11]))
    elif (a_fname == '2G1Ysq_Simul'):
        return ("2G1Ysq", (a_ps[a_idx*3], a_ps[6], a_ps[a_idx*3+1], a_ps[7], a_ps[a_idx*3+2], a_ps[8], a_ps[9]))
    elif (a_fname == '3G1Ysq_Simul'):
        return ("3G1Ysq", (a_ps[a_idx*4], a_ps[8], a_ps[a_idx*4+1], a_ps[9], a_ps[a_idx*4+2], a_ps[10],
                           a_ps[a_idx*4+3], a_ps[11], a_ps[12]))
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_fname)
        return None
