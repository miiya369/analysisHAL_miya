# -*- coding: utf-8 -*-

"""The module to integrate for the folding potential."""

from numpy           import array, inf
from scipy.integrate import quad

def solve_Fint_1D(a_r, a_Integrand_1D, a_Func_1D_opt, a_Params_pot, a_Params_dens):
    """
    The function to integrate (one-dimension) for the folding potential.
    
    For arguments,
    - a_Params_pot [#.param for potential] (1-dim array)
    - a_Params_dens[#.param for density  ] (1-dim array)
    
    return (result of integration, estimated error of integration)
    
    Note1: a_Func_1D_opt  is supposed the pointer of the function for F2(|r|).
    Note2: a_Integrand_1D is supposed the pointer of the function [2Pi Int d_cosP |rP|^2 F1(r-rP) F2(rP)].
    Note3: Integration is defined as f(|r|) = Int_0^inf d_rP Integrand(|r|, |rP| : rho, V).
    """
    
    Params_int = array((a_r, a_Func_1D_opt, a_Params_pot, a_Params_dens))
    
    return quad(a_Integrand_1D, 0.0, inf, args = Params_int)[0]
