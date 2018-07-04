# -*- coding: utf-8 -*-

"""
The module to define the matrix elements for Gauss expansion method 
with folding-potential (Woods-Saxon type).
"""

from numpy import array, empty
from math  import pi, sqrt, exp, erf

from scipy.integrate      import quad
from numpy                import array, inf
from fitting.fitfunc_type import set_fitfunc_from_fname

hbar_c = 197.327053

def set_Ham_mat_gauss_fpot_WS_ij(a_range_i, a_range_j, a_params, a_mass, a_args_dens):
    """
    The function to define a element of Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of multi-ranges Gaussian folding-potential (Woods-Saxon type) ~~~
    
    For arguments,
    - a_range [#.base ] (1-dim ndarray)
    - a_params[#.param] (1-dim ndarray)
    - a_args_dens[3]    (1-dim ndarray)
    
    return: Ham[i,j]
    
    Note1: #.param is got from a_params.
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    """
    l_Nparam = len(a_params)
    
    rho0 = a_args_dens[0]
    RA   = a_args_dens[1]
    aA   = a_args_dens[2]
    
    eRaA = exp(RA/aA)
        
    max_r = 50 # <- You may change here !
    
    hbar_c2_mass = hbar_c**2 / a_mass
        
    nu_ij1 = a_range_i    * a_range_j
    nu_ij2 = a_range_i**2 + a_range_j**2
    tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
    tmp_p = 0.0
    
    for ip in range(0, l_Nparam, 2):
        pk     = a_params[ip + 0]
        qk     = a_params[ip + 1]
        nuijqn = nu_ij1**2 + nu_ij2 * qk**2
        integr = quad(lambda r, Ag: exp(-Ag[0]*r**2) * r**2 / (Ag[1] + exp(r/Ag[2])),
                      0, max_r, args = array((nu_ij2/nuijqn, eRaA, aA)))[0]
        coef   = 8.0*sqrt(2.0)*pi * eRaA * pk*qk**3 * rho0 * pow(nu_ij1 / nuijqn, 1.5)
        tmp_p += coef * integr
    
    return 3.0 * hbar_c2_mass * tmp_d1 / nu_ij2 + tmp_p

from folding_potential.integrand      import integrand_1D_PotOpt_DensWS
from folding_potential.solve_integral import solve_Fint_1D
def set_Ham_mat_fpot_WS_ij(a_range_i, a_range_j, a_params, func_name, a_mass, a_args_dens):
    """
    The function to define a element of Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of optional folding-potential (Woods-Saxon type) ~~~
    
    For arguments,
    - a_params[#.param] (1-dim ndarray)
    - a_args_dens[3]    (1-dim ndarray)
    
    return: Ham[i,j]

    Note: Suppose that a_args_dens[0] = rho_0
    .     Suppose that a_args_dens[1] = R_A
    .     Suppose that a_args_dens[2] = a_A
    """
    rho0 = a_Args_dens[0]
    RA   = a_Args_dens[1]
    aA   = a_Args_dens[2]
    
    l_args = array((integrand_1D_PotOpt_DensWS, set_fitfunc_from_fname(func_name), a_params, array((RA, aA))))
    
    max_r = inf # <- You may change here !
    
    hbar_c2_mass = hbar_c**2 / a_mass
    
    nu_ij1 = a_range_i    * a_range_j
    nu_ij2 = a_range_i**2 + a_range_j**2
    tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
    
    tmp_p  = quad(lambda r, Ag: (solve_Fint_1D(r,Ag[0][0],Ag[0][1],Ag[0][2],Ag[0][3]) * r**2 * exp(-Ag[1]*r**2)), 
                  0, max_r, args = array((l_args, nu_ij2 / nu_ij1**2)))[0]
    tmp_p *= rho0 * 8.0*sqrt(2.0/pi) / pow(nu_ij1, 1.5)
    
    return 3.0 * hbar_c2_mass * tmp_d1 / nu_ij2 + tmp_p

def set_Ham_mat_gauss_fpot_WS(a_range, a_params, a_mass, a_args_dens):
    """
    The function to define a element of Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of multi-ranges Gaussian folding-potential (Woods-Saxon type) ~~~
    
    For arguments,
    - a_range [#.base ] (1-dim ndarray)
    - a_params[#.param] (1-dim ndarray)
    - a_args_dens[3]    (1-dim ndarray)
    
    return: Ham[#.base, #.base] (2-dim ndarray)
    
    Note1: #.base is got from a_range.
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    """
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            r_Ham[i,j] = set_Ham_mat_gauss_fpot_WS_ij(a_range[i], a_range[j], a_params, a_mass, a_Args_dens)
            r_Ham[j,i] = r_Ham[i,j]
    
    return r_Ham

def set_Ham_mat_fpot_WS(a_range, a_params, func_name, a_mass, a_args_dens):
    """
    The function to define a element of Hamiltinian-matrix for Gauss expansion method.
    ~~~ In the case of optional folding-potential (Woods-Saxon type) ~~~
    
    For arguments,
    - a_range [#.base ] (1-dim ndarray)
    - a_params[#.param] (1-dim ndarray)
    - a_args_dens[3]    (1-dim ndarray)
    
    return: Ham[#.base, #.base] (2-dim ndarray)
    
    Note1: #.base is got from a_range.
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    """
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            r_Ham[i,j] = set_Ham_mat_fpot_WS_ij(a_range[i], a_range[j], a_params, func_name, a_mass, a_Args_dens)
            r_Ham[j,i] = r_Ham[i,j]
    
    return r_Ham

from fitting.fitfunc_form import fitfunc_Coulomb
def add_Coulomb_fpot_WS_1(a_range, a_Ham, a_Np, a_args_dens):
    """
    The function to define the Hamiltinian-matrix for Gauss expansion method.
    ~~~ For Coulomb potential of Folding potential with WS-type density ~~~
    
    For arguments,
    - a_range[#.base]         (1-dim ndarray)
    - a_Ham  [#.base, #.base] (2-dim ndarray)
    - a_args_dens[4]          (1-dim ndarray)
    
    return: Ham_C[#.base, #.base] (2-dim ndarray)
    
    Note1: #.base is got from a_range.
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    .      Suppose that a_args_dens[3] = N_A
    """
    if (a_Np == 0):
        return a_Ham
    
    rho0 = a_args_dens[0]
    RA   = a_args_dens[1]
    aA   = a_args_dens[2]
    N_A  = a_args_dens[3]
    
    l_args = array((integrand_1D_PotOpt_DensWS, fitfunc_Coulomb, array((1.0, 1.0)), array((RA, aA))))
    
    max_r = inf # <- You may change here !
    
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1 = a_range[i]    * a_range[j]
            nu_ij2 = a_range[i]**2 + a_range[j]**2
            
            tmp_p  = quad(lambda r, Ag: (solve_Fint_1D(r,Ag[0][0],Ag[0][1],Ag[0][2],Ag[0][3]) * r**2 * exp(-Ag[1]*r**2)),
                          0, max_r, args = array((l_args, nu_ij2 / nu_ij1**2)))[0]
            tmp_p *= rho0 * a_Np * 8.0*sqrt(2.0/pi) / (pow(nu_ij1, 1.5) * N_A)
            
            r_Ham[i,j] = a_Ham[i,j] + tmp_p
            r_Ham[j,i] = r_Ham[i,j]
    
    return r_Ham

from folding_potential.potential_func import fpot_Coulomb_dens_WS
def add_Coulomb_fpot_WS_2(a_range, a_Ham, a_Np, a_args_dens):
    """
    The function to define the Hamiltinian-matrix for Gauss expansion method.
    ~~~ For Coulomb potential of Folding potential with WS-type density ~~~
    
    For arguments,
    - a_range[#.base]         (1-dim ndarray)
    - a_Ham  [#.base, #.base] (2-dim ndarray)
    - a_args_dens[4]          (1-dim ndarray)
    
    return: Ham_C[#.base, #.base] (2-dim ndarray)
    
    Note1: #.base is got from a_range.
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    .      Suppose that a_args_dens[3] = N_A
    """
    if (a_Np == 0):
        return a_Ham
    
    rho0 = a_args_dens[0]
    RA   = a_args_dens[1]
    aA   = a_args_dens[2]
    N_A  = a_args_dens[3]
    
    l_args = array((rho0, RA, aA, a_Np, N_A))
    
    max_r = inf # <- You may change here !
    
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1 = a_range[i]    * a_range[j]
            nu_ij2 = a_range[i]**2 + a_range[j]**2
            
            tmp_p  = quad(lambda r, Ag: (fpot_Coulomb_dens_WS(r, *Ag[0]) * r**2 * exp(-Ag[1]*r**2)), 
                          0, max_r, args = array((l_args, nu_ij2 / nu_ij1**2)))[0]
            tmp_p *= 8.0*sqrt(2.0/pi) / pow(nu_ij1, 1.5)
            
            r_Ham[i,j] = a_Ham[i,j] + tmp_p
            r_Ham[j,i] = r_Ham[i,j]
    
    return r_Ham

def add_Coulomb_fpot_WS_3(a_range, a_Ham, a_Np, a_args_dens):
    """
    The function to define the Hamiltinian-matrix for Gauss expansion method.
    ~~~ For Coulomb potential of Folding potential with WS-type density ~~~
    
    For arguments,
    - a_range[#.base]         (1-dim ndarray)
    - a_Ham  [#.base, #.base] (2-dim ndarray)
    - a_Args_dens[4]          (1-dim ndarray)
    
    return: Ham_C[#.base, #.base] (2-dim ndarray)
    
    Note1: #.base is got from a_range.
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    .      Suppose that a_args_dens[3] = N_A
    """
    if (a_Np == 0):
        return a_Ham
    
    rho0 = a_args_dens[0]
    RA   = a_args_dens[1]
    aA   = a_args_dens[2]
    N_A  = a_args_dens[3]
    
    alp  = 1.43996507808 # = 197.327053 / 137.035999 (= hbar c * alpha)
    CCoe = 8*sqrt(2)*pi * rho0*alp*a_Np/N_A * exp(RA/aA)
    
    max_r = inf # <- You may change here !
    
    l_Nbase = len(a_range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1 = a_range[i]    * a_range[j]
            nu_ij2 = a_range[i]**2 + a_range[j]**2
            
            tmp_p  = quad(lambda r, Ag: (r * erf(Ag[0]*r) * exp(-r/Ag[2]) / (1 + exp((Ag[1] - r) / Ag[2]))), 
                          0, max_r, args = array((sqrt(nu_ij2)/nu_ij1, RA, aA)))[0]
            tmp_p *= (CCoe * pow(nu_ij1/nu_ij2, 1.5))
            
            r_Ham[i,j] = a_Ham[i,j] + tmp_p
            r_Ham[j,i] = r_Ham[i,j]
    
    return r_Ham

def add_Coulomb_fpot_WS(a_range, a_Ham, a_Np, a_args_dens):
    return add_Coulomb_fpot_WS_3(a_range, a_Ham, a_Np, a_args_dens)
