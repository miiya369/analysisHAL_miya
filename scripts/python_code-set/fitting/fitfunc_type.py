# -*- coding: utf-8 -*-

"""The module to verify the fit functions type."""

from fitting.fitfunc_form import *

def get_Nparam_from_fname(a_Fname):
    """
    The function to get the #.param of various fit functions.
    (Fit Fname ==> #.param)
    
    return: #.param
    
    Note: The definition of the Fname is noted in ./definition_function_type.txt
    """
    
    if   (a_Fname == '1P'):
        return 1
    elif (a_Fname == '2P' or a_Fname == '1E' or a_Fname == '1G' or a_Fname == '1CH' or 
          a_Fname == 'SW' or a_Fname == 'Coulomb'):
        return 2
    elif (a_Fname == '3P' or a_Fname == '1SG' or a_Fname == '1Y' or a_Fname == '1Ysq' or
          a_Fname == '1Ytns'):
        return 3
    elif (a_Fname == '4P' or a_Fname == '2E' or a_Fname == '2G' or a_Fname == '2CH'):
        return 4
    elif (a_Fname == '5P'):
        return 5
    elif (a_Fname == '6P'  or a_Fname == '3E' or a_Fname == '3G'   or a_Fname == '3CH' or 
          a_Fname == '2SG' or a_Fname == '2Y' or a_Fname == '2Ysq' or a_Fname == '2Ytns'):
        return 6
    elif (a_Fname == '7P' or a_Fname == '2G1Y' or a_Fname == '2G1Ysq'):
        return 7
    elif (a_Fname == '8P' or a_Fname == '4E' or a_Fname == '4G' or a_Fname == '4CH'):
        return 8
    elif (a_Fname == '9P' or a_Fname == '3SG' or a_Fname == '3Y' or  a_Fname == '3G1Ysq'):
        return 9
    elif (a_Fname == '5E' or a_Fname == '5G' or a_Fname == '5CH' or a_Fname == '2G2Y' or
          a_Fname == '2G1Y1Ysq'):
        return 10
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_Fname)
        return None

def set_fitfunc_from_fname(a_Fname):
    """
    The function to get the fit function.
    (Fit Fname ==> Fit function)
    
    return: Fit function
    
    Note: The definition of the Fname is noted in ./definition_function_type.txt
    """
    
    if   (a_Fname == 'SW'):
        return fitfunc_SW
    elif (a_Fname == 'Coulomb'):
        return fitfunc_Coulomb
    elif (a_Fname == '1E'):
        return fitfunc_1Exp
    elif (a_Fname == '1G'):
        return fitfunc_1G
    elif (a_Fname == '2G'):
        return fitfunc_2G
    elif (a_Fname == '3G'):
        return fitfunc_3G
    elif (a_Fname == '4G'):
        return fitfunc_4G
    elif (a_Fname == '1SG'):
        return fitfunc_1SG
    elif (a_Fname == '2SG'):
        return fitfunc_2SG
    elif (a_Fname == '1Y'):
        return fitfunc_1Y
    elif (a_Fname == '2Y'):
        return fitfunc_2Y
    elif (a_Fname == '1Ysq'):
        return fitfunc_1Ysq
    elif (a_Fname == '2Ysq'):
        return fitfunc_2Ysq
    elif (a_Fname == '1Ytns'):
        return fitfunc_1Ytns
    elif (a_Fname == '2Ytns'):
        return fitfunc_2Ytns
    elif (a_Fname == '2G1Y'):
        return fitfunc_2G1Y
    elif (a_Fname == '2G1Ysq'):
        return fitfunc_2G1Ysq
    elif (a_Fname == '3G1Ysq'):
        return fitfunc_3G1Ysq
    elif (a_Fname == '2G2Y'):
        return fitfunc_2G2Y
    elif (a_Fname == '2G1Y1Ysq'):
        return fitfunc_2G1Y1Ysq
    else:
        print("\nERROR: Invalid (or Have not implemented) fit function type: %s\n" % a_Fname)
        return None
