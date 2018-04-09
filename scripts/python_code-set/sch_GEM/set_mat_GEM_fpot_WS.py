# -*- coding: utf-8 -*-

"""The module to define the matrix elements for Gauss expansion method."""

from numpy import empty, pi, sqrt, exp
from math  import erf
import re

hbar_c = 197.327053

def set_mat_fpot_WS(a_Range, a_Params, a_mu, a_fname, a_Np):
    """
    The function to define the matrix elements.
    
    For arguments,
    - a_Range [#.base ] (1-dim ndarray)
    - a_Params[#.param] (1-dim ndarray)
    
    return: (Ham, Phi)
    """
    if   (re.match('^[1-9]G$', a_fname) is not None):
        r_Ham, r_Phi = set_mat_gauss_fpot_WS(a_Range, a_Params, a_mu[0], a_mu[1:])
    else:
        r_Ham, r_Phi = set_mat_opt_fpot_WS  (a_Range, a_Params, a_fname, a_mu[0], a_mu[1:])
    
    if (a_Np != 0): # For Coulomb potential of Folding potential with WS-type density
        r_Ham = add_Coulomb_fpot_WS(a_Range, r_Ham, a_Np, a_mu[1:])
    
    return (r_Ham, r_Phi)

def set_mat_gauss_fpot_WS_ij(a_range_i, a_range_j, a_Params, a_mu, a_Args_dens):
    """
    The function to define the matrix elements for the folding potential (Woods-Saxon type).
    ~~~ In the case of multi-ranges Gaussian potential ~~~
    
    For arguments,
    - a_Range [#.base ] (1-dim array)
    - a_Params[#.param] (1-dim array)
    - a_Args_dens[3]    (1-dim array)
    
    return: (Ham[i, j], Phi[i, j])
    
    Note1: #.param is got from a_Params.
    Note2: Suppose that a_Args_dens[0] = rho_0
    .      Suppose that a_Args_dens[1] = R_A
    .      Suppose that a_Args_dens[2] = a_A
    """
    
    l_Nparam = len(a_Params)
    
    rho0 = a_Args_dens[0]
    RA   = a_Args_dens[1]
    aA   = a_Args_dens[2]
    
    eRaA = exp(RA/aA)
        
    max_r = 50 # <- You may change here !
    
    hbar_c2_mu = hbar_c**2 / a_mu
        
    nu_ij1 = a_range_i    * a_range_j
    nu_ij2 = a_range_i**2 + a_range_j**2
    tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
    tmp_p = 0.0
    
    for ip in range(0, l_Nparam, 2):
        pk     = a_Params[ip + 0]
        qk     = a_Params[ip + 1]
        nuijqn = nu_ij1**2 + nu_ij2 * qk**2
        integr = quad(lambda r, Ag: exp(-Ag[0]*r**2) * r**2 / (Ag[1] + exp(r/Ag[2])),
                      0, max_r, args = array((nu_ij2/nuijqn, eRaA, aA)))[0]
        coef   = 8.0*sqrt(2.0)*pi * eRaA * pk*qk**3 * rho0 * pow(nu_ij1 / nuijqn, 1.5)
        tmp_p += coef * integr
        
    r_Phi = tmp_d1
    r_Ham = 3.0 * hbar_c2_mu * tmp_d1 / nu_ij2 + tmp_p
    
    return (r_Ham, r_Phi)

from FoldingPotential.Integrand     import integrand_1D_PotOpt_DensWS
from FoldingPotential.SolveIntegral import solve_Fint_1D
def set_mat_opt_fpot_WS_ij(a_range_i, a_range_j, a_Params, a_fname, a_mu, a_Args_dens):
    """
    The function to define the matrix elements for the folding potential (Woods-Saxon type).
    ~~~ In the case of optional potentials ~~~
    
    For arguments,
    - a_Params[#.param] (1-dim array)
    - a_Args_dens[3]    (1-dim array)
    
    return: (Ham[i, j], Phi[i, j])

    Note: Suppose that a_Args_dens[0] = rho_0
    .     Suppose that a_Args_dens[1] = R_A
    .     Suppose that a_Args_dens[2] = a_A
    """
    
    rho0 = a_Args_dens[0]
    RA   = a_Args_dens[1]
    aA   = a_Args_dens[2]
    
    l_Args = array((integrand_1D_PotOpt_DensWS, set_fitfunc_from_fname(a_fname), 
                    a_Params, array((RA, aA))))
    
    l_Func = lambda r, Ag: (solve_Fint_1D(r, Ag[0][0], Ag[0][1], Ag[0][2], Ag[0][3]) 
                            * r**2 * exp(-Ag[1]*r**2))
    
    max_r = inf # <- You may change here !
    
    hbar_c2_mu = hbar_c**2 / a_mu
    
    nu_ij1 = a_range_i    * a_range_j
    nu_ij2 = a_range_i**2 + a_range_j**2
    tmp_d1 = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
    
    tmp_p  = quad(l_Func, 0, max_r, args = array((l_Args, nu_ij2 / nu_ij1**2)))[0]
    tmp_p *= rho0 * 8.0*sqrt(2.0/pi) / pow(nu_ij1, 1.5)
    
    r_Phi = tmp_d1
    r_Ham = 3.0 * hbar_c2_mu * tmp_d1 / nu_ij2 + tmp_p
    
    return (r_Ham, r_Phi)

def set_mat_gauss_fpot_WS(a_Range, a_Params, a_mu, a_Args_dens):
    """
    The function to define the matrix elements for the folding potential (Woods-Saxon type).
    ~~~ In the case of multi-ranges Gaussian potential ~~~
    
    For arguments,
    - a_Range [#.base ] (1-dim array)
    - a_Params[#.param] (1-dim array)
    - a_Args_dens[3]    (1-dim array)
    
    return: (Ham, Phi)
    
    Note1: #.base is got from a_Range.
    Note2: Suppose that a_Args_dens[0] = rho_0
    .      Suppose that a_Args_dens[1] = R_A
    .      Suppose that a_Args_dens[2] = a_A
    """
    
    l_Nbase = len(a_Range)
    
    r_Ham = empty((l_Nbase, l_Nbase))
    r_Phi = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            r_Ham[i, j], r_Phi[i, j] = set_mat_gauss_fpot_WS_ij(a_Range[i], a_Range[j], 
                                                                a_Params, a_mu, a_Args_dens)
            r_Phi[j, i] = r_Phi[i, j]
            r_Ham[j, i] = r_Ham[i, j]
    
    return (r_Ham, r_Phi)

def set_mat_opt_fpot_WS(a_Range, a_Params, a_fname, a_mu, a_Args_dens):
    """
    The function to define the matrix elements for the folding potential (Woods-Saxon type).
    ~~~ In the case of optional potentials ~~~
    
    For arguments,
    - a_Range [#.base ] (1-dim array)
    - a_Params[#.param] (1-dim array)
    - a_Args_dens[3]    (1-dim array)
    
    return: (Ham, Phi)
    
    Note1: #.base is got from a_Range.
    Note2: Suppose that a_Args_dens[0] = rho_0
    .      Suppose that a_Args_dens[1] = R_A
    .      Suppose that a_Args_dens[2] = a_A
    """
    
    l_Nbase = len(a_Range)
    
    r_Ham = empty((l_Nbase, l_Nbase))
    r_Phi = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            r_Ham[i, j], r_Phi[i, j] = set_mat_opt_fpot_WS_ij(a_Range[i], a_Range[j], 
                                                              a_Params, a_fname, a_mu, a_Args_dens)
            r_Phi[j, i] = r_Phi[i, j]
            r_Ham[j, i] = r_Ham[i, j]
    
    return (r_Ham, r_Phi)

from Fitting.FitFunctionForm import fitfunc_Coulomb
def add_Coulomb_fpot_WS_1(a_Range, a_Ham, a_Np, a_Args_dens):
    """
    The function to define the matrix elements for Gauss expansion method.
    ~~~ For Coulomb potential of Folding potential with WS-type density ~~~
    
    For arguments,
    - a_Range[#.base]         (1-dim array)
    - a_Ham  [#.base, #.base] (2-dim array)
    - a_Args_dens[4]          (1-dim array)
    
    return: (Ham)
    
    Note1: #.base is got from a_Range.
    Note2: Suppose that a_Args_dens[0] = rho_0
    .      Suppose that a_Args_dens[1] = R_A
    .      Suppose that a_Args_dens[2] = a_A
    .      Suppose that a_Args_dens[3] = N_A
    """
    
    if (a_Np == 0):
        return a_Ham
    
    rho0 = a_Args_dens[0]
    RA   = a_Args_dens[1]
    aA   = a_Args_dens[2]
    N_A  = a_Args_dens[3]
    
    l_Args = array((integrand_1D_PotOpt_DensWS, fitfunc_Coulomb, 
                    array((1.0, 1.0)), array((RA, aA))))
    
    l_Func = lambda r, Ag: (solve_Fint_1D(r, Ag[0][0], Ag[0][1], Ag[0][2], Ag[0][3])
                            * r**2 * exp(-Ag[1]*r**2))
    
    max_r = inf # <- You may change here !
    
    l_Nbase = len(a_Range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1 = a_Range[i]    * a_Range[j]
            nu_ij2 = a_Range[i]**2 + a_Range[j]**2
            
            tmp_p  = quad(l_Func, 0, max_r, args = array((l_Args, nu_ij2 / nu_ij1**2)))[0]
            tmp_p *= rho0 * a_Np * 8.0*sqrt(2.0/pi) / (pow(nu_ij1, 1.5) * N_A)
            
            r_Ham[i, j] = a_Ham[i, j] + tmp_p
            r_Ham[j, i] = r_Ham[i, j]
    
    return r_Ham

from FoldingPotential.PotentialFunc import fpot_Coulomb_dens_WS
def add_Coulomb_fpot_WS_2(a_Range, a_Ham, a_Np, a_Args_dens):
    """
    The function to define the matrix elements for Gauss expansion method.
    ~~~ For Coulomb potential of Folding potential with WS-type density ~~~
    
    For arguments,
    - a_Range[#.base]         (1-dim array)
    - a_Ham  [#.base, #.base] (2-dim array)
    - a_Args_dens[4]          (1-dim array)
    
    return: (Ham)
    
    Note1: #.base is got from a_Range.
    Note2: Suppose that a_Args_dens[0] = rho_0
    .      Suppose that a_Args_dens[1] = R_A
    .      Suppose that a_Args_dens[2] = a_A
    .      Suppose that a_Args_dens[3] = N_A
    """
    
    if (a_Np == 0):
        return a_Ham
    
    rho0 = a_Args_dens[0]
    RA   = a_Args_dens[1]
    aA   = a_Args_dens[2]
    N_A  = a_Args_dens[3]
    
    l_Args = array((rho0, RA, aA, a_Np, N_A))
    
    l_Func = lambda r, Ag: (fpot_Coulomb_dens_WS(r, *Ag[0]) * r**2 * exp(-Ag[1]*r**2))
    
    max_r = inf # <- You may change here !
    
    l_Nbase = len(a_Range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1 = a_Range[i]    * a_Range[j]
            nu_ij2 = a_Range[i]**2 + a_Range[j]**2
            
            tmp_p  = quad(l_Func, 0, max_r, args = array((l_Args, nu_ij2 / nu_ij1**2)))[0]
            tmp_p *= 8.0*sqrt(2.0/pi) / pow(nu_ij1, 1.5)
            
            r_Ham[i, j] = a_Ham[i, j] + tmp_p
            r_Ham[j, i] = r_Ham[i, j]
    
    return r_Ham

def add_Coulomb_fpot_WS_3(a_Range, a_Ham, a_Np, a_Args_dens):
    """
    The function to define the matrix elements for Gauss expansion method.
    ~~~ For Coulomb potential of Folding potential with WS-type density ~~~
    
    For arguments,
    - a_Range[#.base]         (1-dim array)
    - a_Ham  [#.base, #.base] (2-dim array)
    - a_Args_dens[4]          (1-dim array)
    
    return: (Ham)
    
    Note1: #.base is got from a_Range.
    Note2: Suppose that a_Args_dens[0] = rho_0
    .      Suppose that a_Args_dens[1] = R_A
    .      Suppose that a_Args_dens[2] = a_A
    .      Suppose that a_Args_dens[3] = N_A
    """
    
    if (a_Np == 0):
        return a_Ham
    
    rho0 = a_Args_dens[0]
    RA   = a_Args_dens[1]
    aA   = a_Args_dens[2]
    N_A  = a_Args_dens[3]
    
    alp  = 1.43996507808 # = 197.327053 / 137.035999 (= hbar c * alpha)
    
    l_Func = lambda r, Ag: (r * erf(Ag[0]*r) * exp(-r/Ag[2]) / (1 + exp((Ag[1] - r) / Ag[2])))
    
    Coe = 8*sqrt(2)*pi * rho0*alp*a_Np/N_A * exp(RA/aA)
    
    max_r = inf # <- You may change here !
    
    l_Nbase = len(a_Range)
    r_Ham   = empty((l_Nbase, l_Nbase))
    
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1 = a_Range[i]    * a_Range[j]
            nu_ij2 = a_Range[i]**2 + a_Range[j]**2
            
            tmp_p  = quad(l_Func, 0, max_r, args = array((sqrt(nu_ij2)/nu_ij1, RA, aA)))[0]
            tmp_p *= (Coe * pow(nu_ij1/nu_ij2, 1.5))
            
            r_Ham[i, j] = a_Ham[i, j] + tmp_p
            r_Ham[j, i] = r_Ham[i, j]
    
    return r_Ham

def add_Coulomb_fpot_WS(a_Range, a_Ham, a_Np, a_Args_dens):
    return add_Coulomb_fpot_WS_3(a_Range, a_Ham, a_Np, a_Args_dens)