# -*- coding: utf-8 -*-

"""
The module to solve the Schrodinger equation by solving the differential equation.

~~~ For the case of J+ = 1+ system ~~~
in which S- and D-wave components are mixtured by the tensor force.
"""

hbar_c = 197.327053

from numpy               import array, sqrt, dot, exp, empty, pi
from numpy.linalg        import inv
from Tmatrix.convert_mat import convert_StoT

from misc_QM.special_functions import *

def func_sch_ode_mixLS(a_phi, a_r, a_V_C, a_prm_C, a_V_T, a_prm_T, a_E, a_mu):
    """
    The function for 'odeint' to solve the Schrodinger equation in J+ = 1+ system.
    
    For arguments,
    - a_r, a_E, a_mu   (double scalar)
    - a_phi[4]         (1-dim ndarray) <-- (phi_S, phi_D, d/dr phi_S, d/dr phi_D)
    - a_V_C, a_V_T     (function obj.) <--  central pot, tensor pot
    - a_prm_C[#.param] (1-dim ndarray) <--  parameters of V_C
    - a_prm_T[#.param] (1-dim ndarray) <--  parameters of V_T
    
    Note: The units of potentials, E and mu must be [MeV].
    """
    hbar_c2_mu = hbar_c**2 / a_mu # in unit of [MeV * fm^2]
    
    l_mat = array([[a_V_C(a_r,*a_prm_C) - a_E,
                    a_V_T(a_r,*a_prm_T) * sqrt(8.0)],
                   [a_V_T(a_r,*a_prm_T) * sqrt(8.0),
                    3.0*hbar_c2_mu/a_r**2 + a_V_C(a_r,*a_prm_C) - a_V_T(a_r,*a_prm_T)*2.0 - a_E]
                   ]) * (2.0 / hbar_c2_mu) # in unit of [1/fm^2]
    
    return array([*a_phi[2:], *dot(l_mat, a_phi[:2])])

def solve_sch_diff_mixLS(a_V_C, a_prm_C, a_V_T, a_prm_T, a_mu, a_E, a_rmax):
    """
    The function to solve the Schrodinger equation by solving the differential equation,
    and return the T-matrix.
    
    For arguments,
    - a_rmax, a_mu     (double scalar)
    - a_V_C, a_V_T     (function obj.) <-- central pot, tensor pot
    - a_prm_C[#.param] (1-dim ndarray) <-- parameters of V_C
    - a_prm_T[#.param] (1-dim ndarray) <-- parameters of V_T
    - a_E    [#.data]  (1-dim ndarray) <-- ndarray of energy to calculate
    
    return: T-matrix[#.data, 2, 2] (1-dim ndarray, dtype=complex)
    
    Note1: #.data is got from len(a_E)
    Note2: The units of potentials, E and mu must be [MeV].
    """
    from scipy.integrate import odeint
    
    l_Ndata = len(a_E)
    
    l_rini  = 1e-6 # To avoid zero-div
    l_range = array([l_rini, a_rmax])
    
    # initial values of wave function
    l_iphi  = array([[l_rini**1, l_rini**3, 1.0, 3*l_rini**2],  # S-wave
                     [l_rini**1, l_rini**3, 0.0, 3*l_rini**2]]) # D-wave
    
    r_Tmat  = empty((l_Ndata, 2, 2), dtype=complex)
    
    # Calculation for each energy
    for iE in range(l_Ndata):
        E = a_E[iE] if (a_E[iE] != 0.0) else 1e-10 # To avoid zero-div
        k = sqrt(2.0 * a_mu * E) / hbar_c
        
        ret = array([odeint(func_sch_ode_mixLS, l_iphi[i,:], l_range, 
                            args=(a_V_C, a_prm_C, a_V_T, a_prm_T, E, a_mu))[1,:] for i in range(2)])
        
        hmat = array([[     riccati_hn1(0,k*a_rmax     ), 0.0, riccati_hn2(0,k*a_rmax     ), 0.0],
                      [0.0, riccati_hn1(2,k*a_rmax     ), 0.0, riccati_hn2(2,k*a_rmax     )     ],
                      [     riccati_hn1(0,k*a_rmax,True), 0.0, riccati_hn2(0,k*a_rmax,True), 0.0],
                      [0.0, riccati_hn1(2,k*a_rmax,True), 0.0, riccati_hn2(2,k*a_rmax,True)     ]])
        
        pmat = array([ret[:,0], ret[:,1], ret[:,2]/k, ret[:,3]/k])
        Jmat = dot(inv(hmat), pmat)
        Smat = dot(Jmat[:2,:], inv(Jmat[2:,:]))
        
        r_Tmat[iE,:,:] = convert_StoT(Smat)
    
    return r_Tmat
