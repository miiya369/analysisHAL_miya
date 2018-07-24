# -*- coding: utf-8 -*-

"""The module for scatterings in the Idogata-potential."""

from common.misc     import frange
from numpy           import sqrt, array, sin, append
from scipy.special   import spherical_jn, spherical_yn
from scipy.optimize  import brentq

hbar_c = 197.327053

def tan_d_Idogata3D(a_k, a_V0, a_R0, a_Mu, a_l):
    """
    The function to calculate 'tan[delta_l(k)]' for 3D-Idogata-potential.
    
    return: tan[delta_l(k)]
    
    Note: a_k, a_V0, a_R0 and a_Mu shold be in the unit of [MeV], [MeV], [fm] and [MeV].
    """    
    if (a_k == 0):
        a_k = 10e-10
    
    v0   = a_R0 * sqrt(- 2.0 * a_Mu * a_V0 / hbar_c**2)
    rho  = a_R0 * a_k / hbar_c
    rho1 = sqrt(rho**2 + v0**2)
    
    gamma_l = rho1 * spherical_jn(a_l, rho1, True) / spherical_jn(a_l, rho1) + (a_l + 1.0)
    
    nume = (gamma_l - a_l - 1.0) * spherical_jn(a_l, rho) - rho * spherical_jn(a_l, rho, True)
    deno = (gamma_l - a_l - 1.0) * spherical_yn(a_l, rho) - rho * spherical_yn(a_l, rho, True)
    
    return nume / deno

def Eb_Idogata3D_Swave(a_V0, a_R0, a_Mu):
    """
    The function to calculate binding energies for 3D-Idogata-potential.
    
    return: Eb[#.state] (1-dim ndarray)
    
    Note: a_k, a_V0, a_R0 and a_Mu shold be in the unit of [MeV], [MeV], [fm] and [MeV].
    """
    if (a_V0 > 0):
        return array([0.0])
    
    Ret = array([])
    l_V = -2.0 * a_Mu * a_V0 * a_R0**2 / hbar_c**2
    l_f = lambda x: (sin(x)/x)**2 - 1.0/l_V
    if (l_f(10e-10) < 0):
        return array([0.0])
    
    l_sqrtV  = sqrt(l_V)
    div_x    = l_sqrtV / 10000
    #div_x    = l_sqrtV / 100
    plus_flg = True
    for x in frange(10e-10, l_sqrtV, div_x):
        #print(x, l_f(x))
        
        if (plus_flg):
            if (l_f(x) < 0):
                plus_flg = False
                Ret = append(Ret, ((brentq(l_f, x, x-div_x) * hbar_c / a_R0)**2 + 
                                   2.0*a_Mu*a_V0) / (2.0*a_Mu))
        else:
            if (l_f(x) > 0):
                plus_flg = True
    
    return Ret
