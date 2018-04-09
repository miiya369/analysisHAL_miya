# -*- coding: utf-8 -*-

"""The module to define the matrix elements for Gauss expansion method."""

from numpy import empty, pi, sqrt, exp
from math  import erf
import re

hbar_c = 197.327053

def set_mat(a_Range, a_Params, a_mu, a_fname, a_Np):
    """
    The function to define the matrix elements.
    
    For arguments,
    - a_Range [#.base ] (1-dim ndarray)
    - a_Params[#.param] (1-dim ndarray)
    
    return: (Ham, Phi)
    """
    if   (re.match('^SW$'    , a_fname) is not None):
        r_Ham, r_Phi = set_mat_idogata(a_Range, a_Params, a_mu)
    elif (re.match('^[1-9]G$', a_fname) is not None):
        r_Ham, r_Phi = set_mat_gauss  (a_Range, a_Params, a_mu)
    else:
        r_Ham, r_Phi = set_mat_opt    (a_Range, a_Params, a_fname, a_mu)
    
    if (a_Np != 0): # For Coulomb potential
        r_Ham = add_Coulomb(a_Range, r_Ham, a_Np)
    
    return (r_Ham, r_Phi)

def set_mat_idogata(a_Range, a_Params, a_mu):
    """
    The function to define the matrix elements for Gauss expansion method.
    ~~~ In the case of 3-dim square-well potential ~~~
    
    For arguments,
    - a_Range [#.base] (1-dim array)
    - a_Params[2]      (1-dim array)
    
    return: (Ham, Phi)
    
    Note1: V_0 := a_Params[0], R_0 := a_Params[1] for square-well potential
    Note2: #.base is got from a_Range.
    """
    l_Nbase = len(a_Range)
    
    V_0 = a_Params[0]
    R_0 = a_Params[1]
    
    hbar_c2_mu = hbar_c**2 / a_mu
    
    r_Ham = empty((l_Nbase, l_Nbase))
    r_Phi = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(l_Nbase):
            nu_ij1 = a_Range[i]    * a_Range[j]
            nu_ij2 = a_Range[i]**2 + a_Range[j]**2
            tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
            tmp_d2 = nu_ij2 / nu_ij1**2
            
            r_Phi[i, j] =  tmp_d1
            r_Ham[i, j] = (tmp_d1 * (3.0 * hbar_c2_mu / nu_ij2 +
                                     V_0 * erf(sqrt(tmp_d2) * R_0)) -
                           (pow(2.0, 2.5) * V_0 * R_0 * sqrt(nu_ij1) *
                            exp(-tmp_d2 * R_0**2) / (sqrt(pi) * nu_ij2))
                           )
    return (r_Ham, r_Phi)

def set_mat_gauss(a_Range, a_Params, a_mu):
    """
    The function to define the matrix elements for Gauss expansion method.
    ~~~ In the case of multi-ranges Gaussian potential ~~~
    
    For arguments,
    - a_Range [#.base ] (1-dim array)
    - a_Params[#.param] (1-dim array)

    return: (Ham, Phi)
    
    Note: #.base is got from a_Range. #.param is got from a_Params.
    """
    l_Nbase  = len(a_Range)
    l_Nparam = len(a_Params)
    
    hbar_c2_mu = hbar_c**2 / a_mu
    
    r_Ham = empty((l_Nbase, l_Nbase))
    r_Phi = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(l_Nbase):
            nu_ij1 = a_Range[i]    * a_Range[j]
            nu_ij2 = a_Range[i]**2 + a_Range[j]**2
            tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
            tmp_p = 0.0
            
            for ip in range(0, l_Nparam, 2):
                pk     = a_Params[ip + 0]
                q2k    = a_Params[ip + 1]**2
                tmp_p += pk * pow(2.0 * nu_ij1 * q2k / (nu_ij1**2 + nu_ij2 * q2k), 1.5)
            
            r_Phi[i, j] = tmp_d1
            r_Ham[i, j] = 3.0 * hbar_c2_mu * tmp_d1 / nu_ij2 + tmp_p
    
    return (r_Ham, r_Phi)

from scipy.integrate      import quad
from numpy                import array, inf
from fitting.fitfunc_type import set_fitfunc_from_fname
def set_mat_opt_ij(a_range_i, a_range_j, a_Params, a_fname, a_mu):
    """
    The function to define the matrix elements.
    ~~~ In the case of optional potentials ~~~
    
    For arguments,
    - a_Params[#.param] (1-dim array)
    
    return: (Ham[i, j], Phi[i, j])
    """
    
    func_pot = set_fitfunc_from_fname(a_fname)
    
    l_Func = lambda r, Ag: (func_pot(r, *Ag[0]) * r**2 * exp(-Ag[1]*r**2))
    
    max_r = inf # <- You may change here !
    
    hbar_c2_mu = hbar_c**2 / a_mu
    
    nu_ij1 = a_range_i    * a_range_j
    nu_ij2 = a_range_i**2 + a_range_j**2
    tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
    
    tmp_p  = quad(l_Func, 0, max_r, args = array((a_Params, nu_ij2 / nu_ij1**2)))[0]
    tmp_p *= 8.0*sqrt(2.0/pi) / pow(nu_ij1, 1.5)
    
    r_Phi = tmp_d1
    r_Ham = 3.0 * hbar_c2_mu * tmp_d1 / nu_ij2 + tmp_p
    
    return (r_Ham, r_Phi)

def set_mat_opt(a_Range, a_Params, a_fname, a_mu):
    """
    The function to define the matrix elements.
    ~~~ In the case of optional potentials ~~~
    
    For arguments,
    - a_Range [#.base ] (1-dim array)
    - a_Params[#.param] (1-dim array)
    
    return: (Ham, Phi)
    
    Note: #.base is got from a_Range.
    """
    
    l_Nbase = len(a_Range)
    
    r_Ham = empty((l_Nbase, l_Nbase))
    r_Phi = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            r_Ham[i, j], r_Phi[i, j] = set_mat_opt_ij(a_Range[i], a_Range[j], 
                                                      a_Params, a_fname, a_mu)
            r_Phi[j, i] = r_Phi[i, j]
            r_Ham[j, i] = r_Ham[i, j]
    
    return (r_Ham, r_Phi)

def add_Coulomb(a_Range, a_Ham, a_Np):
    """
    The function to define the matrix elements for Gauss expansion method.
    ~~~ For Coulomb potential ~~~
    
    For arguments,
    - a_Range[#.base]         (1-dim array)
    - a_Ham  [#.base, #.base] (2-dim array)
    
    return: (Ham)
    
    Note: #.base is got from a_Range.
    """
    
    if (a_Np == 0):
        return a_Ham
    
    alp  = 1.43996507808 # = 197.327053 / 137.035999 (= hbar c * alpha)
    CCoe = pow(2.0, 2.5) * alp * a_Np / sqrt(pi)
    
    l_Nbase = len(a_Range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(l_Nbase):
            nu_ij1 = a_Range[i]    * a_Range[j]
            nu_ij2 = a_Range[i]**2 + a_Range[j]**2
            
            r_Ham[i, j] = a_Ham[i, j] + CCoe * sqrt(nu_ij1) / nu_ij2
    
    return r_Ham
