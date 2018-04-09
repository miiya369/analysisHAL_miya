# -*- coding: utf-8 -*-

"""The module to solve the NxN Lippmann-Schwinger equation by the Haftel-Tabakin matrix-inversion method."""

hbar_c = 197.327053

from numpy                     import pi, tan, cos
from numpy.polynomial.legendre import leggauss

def solve_sch_int(a_PotFunc, a_Params, a_Mass, a_Pcut):
    """
    The function to solve the NxN Lippmann-Schwinger equation by the Haftel-Tabakin matrix-inversion method,
    and return the T-matrix.
    
    For arguments,
    - a_Rmax                         (double  scalar )
    - a_PotFunc[#.ch, #.ch]          ( 2-dim  ndarray) <-- list of fit function object
    - a_Params [#.ch, #.ch, #.param] ( 3-dim  ndarray)
    - a_Mass   [#.ch, 2]             ( 2-dim  ndarray) <-- [[1st had, 2nd had], ...]
    
    return: T-matrix[#.data, #.ch, #.ch] (3-dim ndarray, dtype=complex)
    
    Note1: #.ch   is got from len(a_Mass[:,0])
    Note2: #.data is got from len(a_Edata)
    Note3: The units of Potentials, E and Mass must to be [MeV].
    """
    points, weights = leggauss(100)
    
    pj   = a_Pcut *                      tan((pi/4.0)*(points+1.0))
    sj   = a_Pcut * (pi/4.0) * weights / cos((pi/4.0)*(points+1.0))**2
    
    

def pot_pi_pj(a_PotFunc, a_Params, a_pi, a_pj):
    """
    The function to calculate V(p_i, p_j) from V(r).
    
    For arguments,
    - a_PotFunc         (function object)
    - a_Params[#.param] ( 1-dim  ndarray)
    
    return: V(p_i, p_j) [in the unit of MeV]
    
    Note: V(r), p_i and p_j are should be in the unit of [MeV].
    """
    
