# -*- coding: utf-8 -*-

"""The module for Definition of the integrand for the folding potential."""

from math import exp, pi, sqrt, log
from scipy.special import spence

def integrand_1D_PotOpt_DensWS(a_r_int, a_Args):
    """
    The function for Definition of the integrand.
    ~~~  For Woods-Saxon-type nucleus density   ~~~
    ~~~ Potential is optional (sphere function) ~~~
    
    => Defined as Int(|r|, |rP|) = 2Pi Int d_cosP rho_WS(r-rP) [|rP|^2 V(rP)]
    
    For arguments,
    - a_Args[4] (1-dim array)
    
    return (The value of integrand)
    
    Note1: a_r                   := a_Args[0]
    .      a_func_potential      := a_Args[1]
    .      a_Params_pot  (array) := a_Args[2]
    .      a_Params_dens (array) := a_Args[3]
    Note2: Suppose that a_Params_dens[0] = R_A
    .      Suppose that a_Params_dens[1] = a_A
    Note3: Suppose that rho  has no normalized factor (which is rho_0)
    """
    
    # For Debug
    #from numpy           import array
    #from scipy.integrate import quad
    #WS_2D = lambda c0, Ag: (exp ((Ag[2] - sqrt(Ag[0]**2 + Ag[1]**2 + 2.0*Ag[0]*Ag[1]*c0)) / Ag[3]) /
    #                        (exp((Ag[2] - sqrt(Ag[0]**2 + Ag[1]**2 + 2.0*Ag[0]*Ag[1]*c0)) / Ag[3]) + 1.0))
    #ret_integrand = quad(WS_2D, -1.0, 1.0, args = array((a_r_int, a_Args[0], a_Args[3][0], a_Args[3][1])))[0]    
    #return ret_integrand * 2.0*pi * a_r_int**2 * a_Args[1](a_r_int, *a_Args[2])
    
    rprP =     a_Args[0] + a_r_int
    rmrP = abs(a_Args[0] - a_r_int)
    
    exp_p = exp((a_Args[3][0] - rprP) / a_Args[3][1])
    exp_m = exp((a_Args[3][0] - rmrP) / a_Args[3][1])
    
    tmp_d = ((- rprP * log(1+exp_p) + rmrP * log(1+exp_m)) +
             (+     spence(1+exp_p) -     spence(1+exp_m)) * a_Args[3][1])
    
    return tmp_d * 2.0*pi * a_Args[3][1] * a_r_int * a_Args[1](a_r_int, *a_Args[2]) / a_Args[0]

def integrand_1D_PotOpt_DensGauss(a_r_int, a_Args):
    """
    The function for Definition of the integrand.
    ~~~    For Gaussian-type nucleus density    ~~~
    ~~~ Potential is optional (sphere function) ~~~
    
    => Defined as Int(|r|, |rP|) = 2Pi Int d_cosP |rP|^2 rho_G(r-rP) V(rP)
    
    For arguments,
    - a_Args[4] (1-dim array)
    
    return (The value of integrand)
    
    Note1: a_r                   := a_Args[0]
    .      a_func_potential      := a_Args[1]
    .      a_Params_pot  (array) := a_Args[2]
    .      a_Params_dens (array) := a_Args[3]
    Note2: Suppose that a_Params_dens[0] = R_A
    Note3: Suppose that rho has no normalized factor (which is rho_0)
    """
    
    RA2 = a_Args[3][0] ** 2
    
    exp_m = exp(-pow(a_Args[0] - a_r_int, 2) / RA2)
    exp_p = exp(-pow(a_Args[0] + a_r_int, 2) / RA2)
    
    return pi * RA2 * a_r_int * a_Args[1](a_r_int, *a_Args[2]) * (exp_m - exp_p) / a_Args[0]

def integrand_1D_PotGauss_DensOpt(a_r_int, a_Args):
    """
    The function for Definition of the integrand.
    ~~~    For multi-ranges Gaussian-type potential   ~~~
    ~~~ Nucleus density is optional (sphere function) ~~~
    
    => Defined as Int(|r|, |rP|) = 2Pi Int d_cosP |rP|^2 V_MG(r-rP) rho(r)
    
    For arguments,
    - a_Args[4] (1-dim array)
    
    return (The value of integrand)
    
    Note1: a_r                   := a_Args[0]
    .      a_func_density        := a_Args[1]
    .      a_Params_pot  (array) := a_Args[2]
    .      a_Params_dens (array) := a_Args[3]
    Note2: #.param_pot is got from a_Args[2] (a_Params_pot)
    """
    
    l_Nparam_pot = len(a_Args[2])
    
    ret_integrand = 0.0
    for iparam in range(0, l_Nparam_pot, 2):
        pn  = a_Args[2][iparam + 0]
        q2n = a_Args[2][iparam + 1] ** 2
        
        exp_m = exp(-pow(a_Args[0] - a_r_int, 2) / q2n)
        exp_p = exp(-pow(a_Args[0] + a_r_int, 2) / q2n)
        
        ret_integrand += pn * q2n * (exp_m - exp_p)
        
    return ret_integrand * pi * a_r_int * a_Args[1](a_r_int, *a_Args[3]) / a_Args[0]
