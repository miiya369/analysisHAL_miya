# -*- coding: utf-8 -*-

"""The module to convert for T-matrix and S-matrix."""

from numpy import identity

def convert_TtoS(a_Tmat):
    """
    The function to convert from T-matrix to S-matrix.
    
    For arguments,
    - a_Tmat[#.channel, #.channel] (2-dim ndarray)
    
    return: r_Smat[#.channel, #.channel] (4-dim ndarray)
    
    Note1: a_Tmat and r_Smat are supposed to complex-type.
    Note2: #.channel is got from a_Tmat[:,0].
    """
    l_Nch = len(a_Tmat[:,0])
    return (identity(l_Nch) + 2.0j * a_Tmat)

def convert_StoT(a_Smat):
    """
    The function to convert from T-matrix to S-matrix.
    
    For arguments,
    - a_Smat[#.channel, #.channel] (2-dim ndarray)
    
    return: r_Tmat[#.channel, #.channel] (4-dim ndarray)
    
    Note1: a_Tmat and r_Smat are supposed to complex-type.
    Note2: #.channel is got from a_Smat[:,0].
    """
    l_Nch = len(a_Smat[:,0])
    return (0.5j * (identity(l_Nch) - a_Smat))
