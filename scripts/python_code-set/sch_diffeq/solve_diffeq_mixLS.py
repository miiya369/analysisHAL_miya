# -*- coding: utf-8 -*-

"""
The module to solve the Schrodinger equation by solving the differential equation.

~~~ For the case of J+ = 1+ system ~~~
in which S- and D-wave components are mixtured by the tensor force.
"""

hbar_c = 197.327053

from numpy        import array, sqrt, dot, exp, empty, pi
from numpy.linalg import inv

import cmath

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

def solve_sch_diff_mixLS(a_V_C, a_prm_C, a_V_T, a_prm_T, a_mu, a_Edata, a_Rmax):
    """
    The function to solve the Schrodinger equation by solving the differential equation,
    and return the T-matrix.
    
    For arguments,
    - a_Rmax, a_mu     (double scalar)
    - a_V_C, a_V_T     (function obj.) <-- central pot, tensor pot
    - a_prm_C[#.param] (1-dim ndarray) <-- parameters of V_C
    - a_prm_T[#.param] (1-dim ndarray) <-- parameters of V_T
    - a_Edata[#.data]  (1-dim ndarray) <-- ndarray of energy to calculate
    
    return: T-matrix[#.data, 2, 2] (1-dim ndarray, dtype=complex)
    
    Note1: #.data is got from len(a_Edata)
    Note2: The units of potentials, E and mu must be [MeV].
    """
    
    from scipy.integrate import odeint
    
    l_Ndata =   len(a_Edata)
    l_range = array([0, a_Rmax])
    l_iPhi  = array([0.0, 0.0, 1.0, 0.0]) # Set initial values of wave function
    
    r_Tmat  = empty((l_Ndata, 2, 2), dtype=complex)
    
    for iE in range(l_Ndata): # Calculation for each energy
        l_E = a_Edata[iE] # in unit [MeV]
        
        # Set E in unit [MeV], and K in unit [1/fm]
        l_Evec = array([l_E - l_Eth[ich] for ich in range(l_Nch)]); l_Evec[l_Evec==0.0] = 1e-10 # To avoid zero-div
        l_Kvec = array([cmath.sqrt(2.0 * l_Mvec[ich] * l_Evec[ich]) / hbar_c for ich in range(l_Nch)])
        
        l_tmpPhi = odeint(func_sch_ode_mixLS, l_iPhi, l_range, 
                          args=(a_V_C, a_prm_C, a_V_T, a_prm_T, l_E, a_mu))[1]
        
        l_TmpPhi1[:, jch] = l_TmpPhi[:l_Nch ]
        l_TmpPhi2[:, jch] = l_TmpPhi[ l_Nch:]
        
        l_sqKM = sqrt(l_Kvec*l_Mvec)
        l_KexM = exp(-1j *(l_Kvec * a_Rmax - a_l * pi/2.0))
        l_KexP = exp(+1j *(l_Kvec * a_Rmax - a_l * pi/2.0))
        
        l_JmK = dot(diag(l_KexM), l_TmpPhi2 + 1j * dot(diag(l_Kvec), l_TmpPhi1))
        l_JpK = dot(diag(l_KexP), l_TmpPhi2 - 1j * dot(diag(l_Kvec), l_TmpPhi1))
        
        l_Smat = dot(diag(1.0/l_sqKM), dot(l_JmK, dot(inv(l_JpK), diag(l_sqKM))))
        
        r_Tmat[iE, :,:] = 0.5j * (identity(l_Nch) - l_Smat)
    
    return r_Tmat
