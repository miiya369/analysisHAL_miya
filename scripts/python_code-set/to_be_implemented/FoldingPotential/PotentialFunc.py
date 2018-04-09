# -*- coding: utf-8 -*-

"""The module for Definition of the potential functions"""
"""           ~ For the folding-potentials ~           """

from numpy import vectorize, array, sum, sin, cos, pi, exp

sph_j0  = lambda x: sin(x) / x if x != 0 else 1.0
Vsph_j0 = vectorize(sph_j0)

def fpot_Coulomb_dens_FBcoe(a_r, a_Params_FBcoe):
    """
    The function of the folding potential
    - Coulomb potential + Fourier-Bessel Coefficients
    
    For arguments,
    - a_Params_FBcoe[R_cut, a_1, a_2, ...]
    
    return: V(r)
    """
    
    Rcut = a_Params_FBcoe[0]
    a_nu = a_Params_FBcoe[1:]
    nu   = array([i for i in range(1, len(a_nu)+1)])
    
    # 1.43996507808 = 197.327053 / 137.035999 (= hbar c * alpha)
    
    if (a_r <= Rcut):
        Coe = 4.0 * 1.43996507808 * Rcut**2 / pi
        return Coe * sum(a_nu/nu**2 * (Vsph_j0(nu*pi*a_r/Rcut) - cos(nu*pi)))
    else:
        Coe = 4.0 * 1.43996507808 * Rcut**3 / (pi * a_r)
        return Coe * sum(a_nu/nu**2 * (Vsph_j0(nu*pi         ) - cos(nu*pi)))

from mpmath import polylog
def fpot_Coulomb_dens_WS(a_r, rho0, R_A, a_A, N_Z, N_A):
    """
    The function of the folding potential
    - Coulomb potential + WS-type density
    
    return: V(r)
    """
    
    # 1.43996507808 = 197.327053 / 137.035999 (= hbar c * alpha)
    
    Ex          = -exp((R_A - a_r) / a_A)
    Coe         = 4.0 * pi * rho0 * N_Z * 1.43996507808 / N_A
    
    Org_Coulomb = N_Z * 1.43996507808 / a_r
    SFP_Coulomb = Coe * (      a_A**2 * polylog(2, Ex) + 
                         2.0 * a_A**3 * polylog(3, Ex) / a_r)
    
    return Org_Coulomb + SFP_Coulomb

    #Ex2   = -exp(-R_A        / a_A) # For Debug
    #Debug = ((Coe * (a_A**2 *        polylog(2, Ex) + 
    #                 a_A**3 * 2.0 * (polylog(3, Ex) - polylog(3, Ex2)) / a_r +
    #                 R_A * (pi**2 * a_A**2 + R_A**2) / (3.0 * a_r))) - 
    #         (Org_Coulomb + SFP_Coulomb))
    #
    #return (Org_Coulomb + SFP_Coulomb, Org_Coulomb, SFP_Coulomb, Debug)
