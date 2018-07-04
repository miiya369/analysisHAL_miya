# -*- coding: utf-8 -*-

"""The module to define the matrix elements for Gauss expansion method."""

from numpy import array, empty
from math  import pi, sqrt, exp, erf

hbar_c = 197.327053

def set_Phi_mat(a_range):
    """
    The function to define the Phi-matrix for Gauss expansion method.
    
    For arguments,
    - a_range[#.base] (1-dim ndarray)
    
    return: Phi[#.base, #.base] (2-dim ndarray)
    
    Note: #.base is got from a_range.
    """
    l_Nbase = len(a_range)
    return array([[pow(2.0*a_range[i]*a_range[j] / (a_range[i]**2+a_range[j]**2), 1.5) 
                   for j in range(l_Nbase)] for i in range(l_Nbase)])

def set_Ham_mat_idogata(a_range, a_params, a_mass):
    """
    The function to define the Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of 3-dim square-well potential ~~~
    
    For arguments,
    - a_range [#.base] (1-dim ndarray)
    - a_params[2]      (1-dim ndarray)
    
    return: Ham[#.base, #.base] (2-dim ndarray)
    
    Note1: V_0 := a_params[0], R_0 := a_params[1] for square-well potential
    Note2: #.base is got from a_range.
    """
    l_Nbase = len(a_range)
    V_0     = a_params[0]
    R_0     = a_params[1]
    
    hbar_c2_mass = hbar_c**2 / a_mass
    
    r_Ham = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(l_Nbase):
            nu_ij1 = a_range[i]    * a_range[j]
            nu_ij2 = a_range[i]**2 + a_range[j]**2
            tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
            tmp_d2 = nu_ij2 / nu_ij1**2
            
            r_Ham[i,j] = (tmp_d1 * (3.0 * hbar_c2_mass / nu_ij2 + V_0 * erf(sqrt(tmp_d2) * R_0)) -
                          (pow(2.0, 2.5) * V_0 * R_0 * sqrt(nu_ij1) * exp(-tmp_d2 * R_0**2) / (sqrt(pi) * nu_ij2))
                          )
    return r_Ham

def set_Ham_mat_gauss(a_range, a_params, a_mass):
    """
    The function to define the Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of multi-ranges Gaussian potential ~~~
    
    For arguments,
    - a_range [#.base ] (1-dim ndarray)
    - a_params[#.param] (1-dim ndarray)
    
    return: Ham[#.base, #.base] (2-dim ndarray)
    
    Note: #.base is got from a_range. #.param is got from a_params.
    """
    l_Nbase  = len(a_range)
    l_Nparam = len(a_params)
    
    hbar_c2_mass = hbar_c**2 / a_mass
    
    r_Ham = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(l_Nbase):
            nu_ij1 = a_range[i]    * a_range[j]
            nu_ij2 = a_range[i]**2 + a_range[j]**2
            tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
            tmp_p = 0.0
            
            for ip in range(0, l_Nparam, 2):
                pk     = a_params[ip + 0]
                q2k    = a_params[ip + 1]**2
                tmp_p += pk * pow(2.0 * nu_ij1 * q2k / (nu_ij1**2 + nu_ij2 * q2k), 1.5)
            
            r_Ham[i,j] = 3.0 * hbar_c2_mass * tmp_d1 / nu_ij2 + tmp_p
    
    return r_Ham

from scipy.integrate      import quad
from numpy                import array, inf
from fitting.fitfunc_type import set_fitfunc_from_fname
def set_Ham_mat_ij(a_range_i, a_range_j, a_params, func_name, a_mass):
    """
    The function to define a element of Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of optional potentials ~~~
    
    For arguments,
    - a_params[#.param] (1-dim ndarray)
    
    return: Ham[i,j]
    """    
    func_pot = set_fitfunc_from_fname(func_name)
    
    max_r = inf # <- You may change here !
    
    hbar_c2_mass = hbar_c**2 / a_mass
    
    nu_ij1 = a_range_i    * a_range_j
    nu_ij2 = a_range_i**2 + a_range_j**2
    tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
    
    tmp_p  = quad(lambda r, Ag: (func_pot(r, *Ag[0]) * r**2 * exp(-Ag[1]*r**2)), 
                  0, max_r, args = array((a_params, nu_ij2 / nu_ij1**2)))[0]
    tmp_p *= 8.0*sqrt(2.0/pi) / pow(nu_ij1, 1.5)
    
    return 3.0 * hbar_c2_mass * tmp_d1 / nu_ij2 + tmp_p

def set_Ham_mat(a_range, a_params, func_name, a_mass):
    """
    The function to define a element of Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of optional potentials ~~~
    
    For arguments,
    - a_range [#.base ] (1-dim ndarray)
    - a_params[#.param] (1-dim ndarray)
    
    return: Ham[#.base, #.base] (2-dim ndarray)
    
    Note: #.base is got from a_range.
    """    
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            r_Ham[i,j] = set_Ham_mat_ij(a_range[i], a_range[j], a_params, func_name, a_mass)
            r_Ham[j,i] = r_Ham[i,j]
    
    return r_Ham

def add_Coulomb(a_range, a_Ham, a_Np):
    """
    The function to define the Hamiltinian-matrix for Gauss expansion method.
    ~~~ For Coulomb potential ~~~
    
    For arguments,
    - a_range[#.base]         (1-dim ndarray)
    - a_Ham  [#.base, #.base] (2-dim ndarray)
    
    return: Ham_C[#.base, #.base] (2-dim ndarray)
    
    Note: #.base is got from a_range.
    """
    if (a_Np == 0):
        return a_Ham
    
    alp  = 1.43996507808 # = 197.327053 / 137.035999 (= hbar c * alpha)
    CCoe = pow(2.0, 2.5) * alp * a_Np / sqrt(pi)
    
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(l_Nbase):
            nu_ij1 = a_range[i]    * a_range[j]
            nu_ij2 = a_range[i]**2 + a_range[j]**2
            
            r_Ham[i,j] = a_Ham[i,j] + CCoe * sqrt(nu_ij1) / nu_ij2
    
    return r_Ham
