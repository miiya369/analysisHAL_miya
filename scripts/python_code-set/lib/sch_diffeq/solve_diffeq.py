# -*- coding: utf-8 -*-

"""The module to solve the NxN Schrodinger equation by solving the differential equation."""

hbar_c = 197.327053

from numpy                     import array, dot, diag, sqrt, exp, empty, pi
from numpy.linalg              import inv
from Tmatrix.convert_mat       import convert_StoT
from misc_QM.special_functions import riccati_hn1, riccati_hn2

import cmath

def func_Helmholtz_ode(a_phi, a_rho, a_l):
    """
    The function for 'odeint' to solve the Helmholtz equation.
    
    For arguments,
    - a_rho    (double scalar)
    - a_phi[2] (1-dim ndarray)
    - a_l      (int    scalar)
    
    Note: The all of variables are assumed as dimension-less
    """
    return array([a_phi[1], (a_l*(a_l+1)/a_rho**2 - 1) * a_phi[0]])

def func_sch_ode(a_phi, a_r, a_V, a_prm, V_dummy, prm_dummy, a_E, a_mu, a_l):
    """
    The function for 'odeint' to solve the NxN Schrodinger equation.
    
    For arguments,
    - a_phi[2 * #.ch]            (1-dim ndarray)
    - a_r                        (double scalar)
    - a_V  [#.ch, #.ch]          (2-dim ndarray) <-- list of fit function object
    - a_prm[#.ch, #.ch, #.param] (3-dim ndarray)
    - a_mu [#.ch]                (1-dim ndarray) <-- a_mu = [mu1, mu2, ...] (Not 2*mu !!)
    - a_E  [#.ch]                (1-dim ndarray) <-- a_E  = [ E0,  E1, ...]
    - a_l  [#.ch]                (1-dim ndarray) <-- a_l  = [ l1,  l2, ...]
    
    Note1 The units of V, E and mu must be [MeV].
    Note2: #.ch is got from a_E.
    """    
    l_Nch = len(a_E)
    
    # in unit [MeV]
    l_Pot = array([[a_V[ich,jch](a_r,*a_prm[ich,jch]) for jch in range(l_Nch)] for ich in range(l_Nch)])
    
    # in unit [1/fm^2]
    l_Amat = (1.0/hbar_c**2) * dot(diag(2.0*a_mu), (l_Pot-diag(a_E))) + diag(a_l*(a_l+1.0)/a_r**2)
    
    return array([*a_phi[l_Nch:], *dot(l_Amat, a_phi[:l_Nch])])

def func_sch_ode_NLO(a_phi, a_r, a_V_LO, a_prm_LO, a_V_NLO, a_prm_NLO, a_E, a_mu, a_l):
    """
    The function for 'odeint' to solve the NxN Schrodinger equation.
    (For LO-NLO potential).
    
    For arguments,
    - a_phi    [2 * #.ch]            (1-dim ndarray)
    - a_r                            (double scalar)
    - a_V_LO   [#.ch, #.ch]          (2-dim ndarray) <-- list of fit function object
    - a_V_NLO  [#.ch, #.ch]          (2-dim ndarray) <-- list of fit function object
    - a_prm_LO [#.ch, #.ch, #.param] (3-dim ndarray)
    - a_prm_NLO[#.ch, #.ch, #.param] (3-dim ndarray)
    - a_mu     [#.ch]                (1-dim ndarray) <-- a_mu = [mu1, mu2, ...] (Not 2*mu !!)
    - a_E      [#.ch]                (1-dim ndarray) <-- a_E  = [ E0,  E1, ...]
    - a_l      [#.ch]                (1-dim ndarray) <-- a_l  = [ l1,  l2, ...]
    
    Note1: The units of V_LO, V_NLO, E and mu must be [MeV].
    Note2: #.ch is got from a_E.
    """
    l_Nch = len(a_E)
    
    # in unit [MeV]
    l_Pot_LO  = array([[a_V_LO [ich,jch](a_r, *a_prm_LO [ich,jch]) for jch in range(l_Nch)] for ich in range(l_Nch)])
    
    # in unit [1/MeV]
    l_Pot_NLO = array([[a_V_NLO[ich,jch](a_r, *a_prm_NLO[ich,jch]) for jch in range(l_Nch)] for ich in range(l_Nch)])
    
    # in unit [1/fm^2]
    l_Amat = (1.0/hbar_c**2) * dot(inv(diag(1.0/(2.0*a_mu))-l_Pot_NLO), (l_Pot_LO-diag(a_E))) + diag(a_l*(a_l+1.0)/a_r**2)
    
    return array([*a_phi[l_Nch:], *dot(l_Amat, a_phi[:l_Nch])])

def solve_sch_diff(a_iphi, a_V, a_prm, a_m, a_E, a_l, a_rini, a_rmax):
    """
    The function to solve the NxN Schrodinger equation by solving the differential equation,
    and return the T-matrix.
    
    For arguments,
    - a_iphi[#.ch, #.ch * 2]      (2-dim ndarray) <-- [[iphi1, iphi2, ..., d/dr iphi1, d/dr iphi2, ...], ...]
    - a_V   [#.ch, #.ch]          (2-dim ndarray) <-- list of fit function object
    - a_prm [#.ch, #.ch, #.param] (3-dim ndarray)
    - a_m   [#.ch, 2]             (2-dim ndarray) <-- [[1st had, 2nd had], ...]
    - a_E   [#.data]              (1-dim ndarray) <-- ndarray of energy to calculate
    - a_l   [#.ch]                (1-dim ndarray) <-- [l1, l2, ...] (integer)
    - a_rini, a_rmax              (double scala)
    
    return: T-matrix[#.data, #.ch, #.ch] (3-dim ndarray, dtype=complex)
    
    Note1: #.ch   is got from len(a_m[:,0])
    Note2: #.data is got from len(a_E)
    Note3: The units of Potentials, E and m must to be [MeV]
    Note4: The channel masses are assumed to arrange in ascending order
    """
    return solve_sch_diff_NLO(a_iphi, a_V, a_prm, None, None, a_m, a_E, a_l, a_rini, a_rmax)

def solve_sch_diff_NLO(a_iphi, a_V_LO, a_prm_LO, a_V_NLO, a_prm_NLO, a_m, a_E, a_l, a_rini, a_rmax):
    """
    The function to solve the NxN Schrodinger equation by solving the differential equation,
    (For LO-NLO potential) and return the T-matrix.
    
    For arguments,
    - a_iphi   [#.ch, #.ch * 2]      (2-dim ndarray) <-- [[iphi1, iphi2, ..., d/dr iphi1, d/dr iphi2, ...], ...]
    - a_V_LO   [#.ch, #.ch]          (2-dim ndarray) <-- list of fit function object
    - a_V_NLO  [#.ch, #.ch]          (2-dim ndarray) <-- list of fit function object
    - a_prm_LO [#.ch, #.ch, #.param] (3-dim ndarray)
    - a_prm_NLO[#.ch, #.ch, #.param] (3-dim ndarray)
    - a_m      [#.ch, 2]             (2-dim ndarray) <-- [[1st had, 2nd had], ...]
    - a_E      [#.data]              (1-dim ndarray) <-- ndarray of energy to calculate
    - a_l      [#.ch]                (1-dim ndarray) <-- [l1, l2, ...] (integer)
    - a_rini, a_rmax                 (double scalar)
    
    return: T-matrix[#.data, #.ch, #.ch] (3-dim ndarray, dtype=complex)
    
    Note1: #.ch   is got from len(a_m[:,0])
    Note2: #.data is got from len(a_E)
    Note3: The units of Potentials, E and m must to be [MeV]
    Note4: The channel masses are assumed to arrange in ascending order
    """
    from scipy.integrate import odeint
    
    l_Nch   = len(a_m[:,0])
    l_Ndata = len(a_E)
    
    # Set vectors for mu and E_th in unit [MeV]
    l_Mvec = array([(a_m[ich,0]*a_m[ich,1]) / (a_m[ich,0]+a_m[ich,1]) for ich in range(l_Nch)])
    l_Eth  = array([(a_m[ich,0]+a_m[ich,1]) - (a_m[  0,0]+a_m[  0,1]) for ich in range(l_Nch)])
    
    l_range  = array([a_rini, a_rmax])
    r_Tmat   = empty((l_Ndata, l_Nch, l_Nch), dtype=complex)
    func_ode = func_sch_ode if (a_V_NLO is a_prm_NLO is None) else func_sch_ode_NLO
    
    # Calculation for each energy
    for iE in range(l_Ndata):
        E = a_E[iE] # in unit [MeV]
        
        # Set vectors for E in unit [MeV], and K in unit [1/fm]
        Evec = array([E - l_Eth[ich] for ich in range(l_Nch)]); Evec[Evec==0.0] = 1e-10 # To avoid zero-div
        Kvec = array([cmath.sqrt(2.0 * l_Mvec[ich] * Evec[ich]) / hbar_c for ich in range(l_Nch)], dtype=complex)
        
        ophi = array([odeint(func_ode, a_iphi[i,:], l_range,
                             args=(a_V_LO, a_prm_LO, a_V_NLO, a_prm_NLO, Evec, l_Mvec, a_l))[1,:]
                      for i in range(l_Nch)]).T
        
        hnL = array([[riccati_hn1(a_l, Kvec*a_rmax      ), riccati_hn2(a_l, Kvec*a_rmax      )],
                     [riccati_hn1(a_l, Kvec*a_rmax, True), riccati_hn2(a_l, Kvec*a_rmax, True)]])
        
        J_K = array([dot(diag(1.0 / (hnL[0,1-i] / hnL[0,i] - hnL[1,1-i] / hnL[1,i])), 
                         (dot(diag(Kvec/hnL[0,i]), ophi[:l_Nch ,:]) - 
                          dot(diag(1.0 /hnL[1,i]), ophi[ l_Nch:,:]))) for i in range(2)])
        
        sqKM = sqrt(Kvec*l_Mvec)
        Smat = dot(diag(1.0/sqKM), dot(J_K[1], dot(inv(J_K[0]), diag(sqKM))))
        
        r_Tmat[iE,:,:] = convert_StoT(Smat)
    
    return r_Tmat
