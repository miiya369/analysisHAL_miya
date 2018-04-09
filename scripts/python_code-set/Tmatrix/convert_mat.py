# -*- coding: utf-8 -*-

"""The module to convert for T-matrix and S-matrix."""

def convert_TtoS(a_Tmat):
    """
    The function to convert from T-matrix to S-matrix.
    
    For arguments,
    - a_Tmat[#.conf, #.data, #.channel, #.channel] (4-dim ndarray)
    
    Note1: a_Tmat is supposed to complex-type.
    Note2: #.channel is got from a_Tmat[0,0,:,0].
    """
    
    a_Tmat *= complex(0.0, 2.0)
    for ich in range(len(a_Tmat[0,0,:,0])):
        a_Tmat[:,:,ich,ich] += complex(1.0, 0.0)

def convert_TtoS_single(a_Tmat):
    """
    The function to convert from T-matrix to S-matrix.
    
    For arguments,
    - a_Tmat[#.conf, #.data] (2-dim ndarray)
    
    Note: a_Tmat is supposed to complex-type.
    """
    
    a_Tmat *= complex(0.0, 2.0)
    a_Tmat += complex(1.0, 0.0)
