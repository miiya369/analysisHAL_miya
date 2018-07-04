# -*- coding: utf-8 -*-

"""The module to input/output text data."""

from numpy import loadtxt, savetxt

def input_text_data(a_ifname, verbose_flg = True):
    """
    The function to input the text data.
    
    return: Data[#.row, #.column]
    
    Note: The text data is assumed to construct as (#.row) x (#.column) float data.
    """
    
    r_Data = loadtxt(a_ifname)
    
    if (verbose_flg):
        print("# Successful to input text data")
        print("# N.row    = %d" % len(r_Data[:,0]))
        print("# N.column = %d" % len(r_Data[0,:]))
    
    return r_Data
    
def output_text_data(a_ofname, a_Data, verbose_flg = True):
    """
    The function to output as the (#.row) x (#.column) float text data.
    
    For arguments,
    - a_Data[#.row, #.column] (2-dim array)
    
    Note: #.row and #.column are got from a_Data.
    """
    
    savetxt(a_ofname, a_Data, fmt='%1.16e '*(len(a_Data[0,:])))
    
    if (verbose_flg):
        print("# Successful to output text data.")
