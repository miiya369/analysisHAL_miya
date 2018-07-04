# -*- coding: utf-8 -*-

"""The module to verify the fit functions type."""

from fitting.fitfunc_form import *

def get_Nparam_from_fname(a_fname):
    """
    The function to get the #.param of various fit functions.
    (Fit fname ==> #.param)
    
    return: #.param
    
    Note: The definition of the fname is noted in ./definition_function_type.txt
    """
    
    if   (a_fname == '1P'):
        return 1
    elif (a_fname == '2P' or a_fname == '1E' or a_fname == '1G' or a_fname == '1CH' or 
          a_fname == 'SW' or a_fname == 'Coulomb'):
        return 2
    elif (a_fname == '3P' or a_fname == '1SG' or a_fname == '1Y' or a_fname == '1Ysq' or
          a_fname == '1Ytns'):
        return 3
    elif (a_fname == '4P' or a_fname == '2E' or a_fname == '2G' or a_fname == '2CH'):
        return 4
    elif (a_fname == '5P' or a_fname == '1G1Y'):
        return 5
    elif (a_fname == '6P'  or a_fname == '3E' or a_fname == '3G'   or a_fname == '3CH' or 
          a_fname == '2SG' or a_fname == '2Y' or a_fname == '2Ysq' or a_fname == '2Ytns'):
        return 6
    elif (a_fname == '7P' or a_fname == '2G1Y' or a_fname == '2G1Ysq'):
        return 7
    elif (a_fname == '8P' or a_fname == '4E' or a_fname == '4G' or a_fname == '4CH'):
        return 8
    elif (a_fname == '9P' or a_fname == '3SG' or a_fname == '3Y' or  a_fname == '3G1Ysq'):
        return 9
    elif (a_fname == '5E' or a_fname == '5G' or a_fname == '5CH' or a_fname == '2G2Y' or
          a_fname == '2G1Y1Ysq'):
        return 10
    elif (a_fname == '12G'):
        return 24
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_fname)
        return None

def set_fitfunc_from_fname(a_fname):
    """
    The function to get the fit function.
    (Fit fname ==> Fit function)
    
    return: Fit function
    
    Note: The definition of the fname is noted in ./definition_function_type.txt
    """
    
    if   (a_fname == 'SW'):
        return fitfunc_SW
    elif (a_fname == 'Coulomb'):
        return fitfunc_Coulomb
    elif (a_fname == '1E'):
        return fitfunc_1Exp
    elif (a_fname == '1G'):
        return fitfunc_1G
    elif (a_fname == '2G'):
        return fitfunc_2G
    elif (a_fname == '3G'):
        return fitfunc_3G
    elif (a_fname == '4G'):
        return fitfunc_4G
    elif (a_fname == '12G'):
        return fitfunc_12G
    elif (a_fname == '1SG'):
        return fitfunc_1SG
    elif (a_fname == '2SG'):
        return fitfunc_2SG
    elif (a_fname == '1Y'):
        return fitfunc_1Y
    elif (a_fname == '2Y'):
        return fitfunc_2Y
    elif (a_fname == '1Ysq'):
        return fitfunc_1Ysq
    elif (a_fname == '2Ysq'):
        return fitfunc_2Ysq
    elif (a_fname == '1Ytns'):
        return fitfunc_1Ytns
    elif (a_fname == '2Ytns'):
        return fitfunc_2Ytns
    elif (a_fname == '1G1Y'):
        return fitfunc_1G1Y
    elif (a_fname == '2G1Y'):
        return fitfunc_2G1Y
    elif (a_fname == '2G1Ysq'):
        return fitfunc_2G1Ysq
    elif (a_fname == '3G1Ysq'):
        return fitfunc_3G1Ysq
    elif (a_fname == '2G2Y'):
        return fitfunc_2G2Y
    elif (a_fname == '2G1Y1Ysq'):
        return fitfunc_2G1Y1Ysq
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_fname)
        return None
