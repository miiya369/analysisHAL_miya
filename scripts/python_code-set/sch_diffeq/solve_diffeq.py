# -*- coding: utf-8 -*-

"""The module to solve the NxN Schrodinger equation by solving the differential equation."""

hbar_c = 197.327053

from numpy        import array, dot, diag, identity, zeros, sqrt, exp, empty, pi
from numpy.linalg import inv

import cmath

def func_sch_ode(a_Phi, a_r, a_PotFunc, a_Params, a_E, a_Mu, a_l = 0):
    """
    The function for 'odeint' to solve the NxN Schrodinger equation.
    
    For arguments,
    - a_r                            (double  scalar )
    - a_Phi    [2 * #.ch]            ( 1-dim  ndarray)
    - a_PotFunc[#.ch, #.ch]          ( 2-dim  ndarray) <-- list of fit function object
    - a_Params [#.ch, #.ch, #.param] ( 3-dim  ndarray)
    - a_Mu     [#.ch]                ( 1-dim  ndarray) <-- a_Mu = [Mu1, Mu2, ...] (Not 2*Mu !!)
    - a_E      [#.ch]                ( 1-dim  ndarray) <-- a_E  = [ E0,  E1, ...]
    
    Note1: The units of Pot, E and Mu must be [MeV].
    Note2: #.ch is got from a_E.
    """
    
    l_Nch  = len(a_E)
    l_Pot  = array([[a_PotFunc[ich, jch](a_r, *a_Params[ich, jch])
                     for jch in range(l_Nch)] for ich in range(l_Nch)])
    
    l_Amat = (1.0/hbar_c**2) * dot(diag(2.0*a_Mu), (l_Pot - diag(a_E))) # in unit [1/fm^2]
    #+ a_l*(a_l+1.0)/a_r**2 * identity(l_Nch)
    
    return array([*a_Phi[l_Nch:], *dot(l_Amat, a_Phi[:l_Nch])])

def func_sch_ode_NLO(a_Phi, a_r, a_PotFunc__LO, a_Params__LO, a_PotFunc_NLO, a_Params_NLO, 
                     a_E, a_Mu, a_l = 0):
    """
    The function for 'odeint' to solve the NxN Schrodinger equation.
    (For LO-NLO potential).
    
    For arguments,
    - a_r                                               (double  scalar )
    - a_Phi                       [2 * #.ch]            ( 1-dim  ndarray)
    - a_PotFunc__LO, a_PotFunc_NLO[#.ch, #.ch]          ( 2-dim  ndarray) <-- list of fit function object
    - a_Params__LO,  a_Params_NLO [#.ch, #.ch, #.param] ( 3-dim  ndarray)
    - a_Mu                        [#.ch]                ( 1-dim  ndarray) <-- a_Mu = [Mu1, Mu2, ...] (Not 2*Mu !!)
    - a_E                         [#.ch]                ( 1-dim  ndarray) <-- a_E  = [ E0,  E1, ...]
    
    Note1: The units of Pot__LO, Pot_NLO, E and Mu must be [MeV].
    Note2: #.ch is got from a_E.
    """
    
    l_Nch     = len(a_E)
    l_Pot__LO = array([[a_PotFunc__LO[ich, jch](a_r, *a_Params__LO[ich, jch])
                        for jch in range(l_Nch)] for ich in range(l_Nch)])     # in unit [  MeV]
    
    l_Pot_NLO = array([[a_PotFunc_NLO[ich, jch](a_r, *a_Params_NLO[ich, jch])
                        for jch in range(l_Nch)] for ich in range(l_Nch)])     # in unit [1/MeV]
        
    l_Amat    = (1.0/hbar_c**2) * dot(inv(diag(1.0/(2.0*a_Mu)) - l_Pot_NLO), 
                                      (l_Pot__LO - diag(a_E)))                 # in unit [1/fm^2]
    #+ a_l*(a_l+1.0)/a_r**2 * identity(l_Nch)
    
    return array([*a_Phi[l_Nch:], *dot(l_Amat, a_Phi[:l_Nch])])

def solve_sch_diff(a_PotFunc, a_Params, a_Mass, a_Edata, a_Rmax, a_l = 0):
    """
    The function to solve the NxN Schrodinger equation by solving the differential equation,
    and return the T-matrix.
    
    For arguments,
    - a_Rmax                         (double  scalar )
    - a_PotFunc[#.ch, #.ch]          ( 2-dim  ndarray) <-- list of fit function object
    - a_Params [#.ch, #.ch, #.param] ( 3-dim  ndarray)
    - a_Mass   [#.ch, 2]             ( 2-dim  ndarray) <-- [[1st had, 2nd had], ...]
    - a_Edata  [#.data]              ( 1-dim  ndarray) <-- ndarray of energy to calculate
    
    return: T-matrix[#.data, #.ch, #.ch] (3-dim ndarray, dtype=complex)
    
    Note1: #.ch   is got from len(a_Mass[:,0])
    Note2: #.data is got from len(a_Edata)
    Note3: The units of Potentials, E and Mass must to be [MeV].
    """
    
    return solve_sch_diff_NLO(a_PotFunc, a_Params, None, None, a_Mass, a_Edata, a_Rmax, a_l)

def solve_sch_diff_NLO(a_PotFunc__LO, a_Params__LO, a_PotFunc_NLO, a_Params_NLO, 
                       a_Mass, a_Edata, a_Rmax, a_l = 0):
    """
    The function to solve the NxN Schrodinger equation by solving the differential equation,
    (For NLO potential) and return the T-matrix.
    
    For arguments,
    - a_Rmax                                            (double  scalar )
    - a_PotFunc__LO, a_PotFunc_NLO[#.ch, #.ch]          ( 2-dim  ndarray) <-- list of fit function object
    - a_Params__LO,  a_Params_NLO [#.ch, #.ch, #.param] ( 3-dim  ndarray)
    - a_Mass                      [#.ch, 2]             ( 2-dim  ndarray) <-- [[1st had, 2nd had], ...]
    - a_Edata                     [#.data]              ( 1-dim  ndarray) <-- ndarray of energy to calculate
    
    return: T-matrix[#.data, #.ch, #.ch] (3-dim ndarray, dtype=complex)
    
    Note1: #.ch   is got from len(a_Mass[:,0])
    Note2: #.data is got from len(a_Edata)
    Note3: The units of Potentials, E and Mass must to be [MeV].
    """
    
    from scipy.integrate import odeint
    
    l_Nch   = len(a_Mass[:,0])
    l_Ndata = len(a_Edata)
    l_Range = array([0, a_Rmax])
    
    # Set Mu and E_th in unit [MeV]
    l_Mvec = array([(a_Mass[ich,0]*a_Mass[ich,1]) / (a_Mass[ich,0]+a_Mass[ich,1]) for ich in range(l_Nch)])
    l_Eth  = array([(a_Mass[ich,0]+a_Mass[ich,1]) - (a_Mass[  0,0]+a_Mass[  0,1]) for ich in range(l_Nch)])
    
    l_InitPhi = zeros((2*l_Nch, l_Nch), dtype=float) # Set initial values of wave function
    for ich in range(l_Nch):
        l_InitPhi[ich+l_Nch, ich] = 1.0
    
    r_Tmat    = empty((l_Ndata, l_Nch, l_Nch), dtype=complex)
    l_TmpPhi1 = empty((l_Nch, l_Nch), dtype=float)
    l_TmpPhi2 = empty((l_Nch, l_Nch), dtype=float)
    for iE in range(l_Ndata): # Calculation for each energy
        l_E = a_Edata[iE] # in unit [MeV]
        
        # Set E in unit [MeV], and K in unit [1/fm]
        l_Evec = array([l_E - l_Eth[ich] for ich in range(l_Nch)]); l_Evec[l_Evec==0.0] = 1e-10 # To avoid zero-div
        l_Kvec = array([cmath.sqrt(2.0 * l_Mvec[ich] * l_Evec[ich]) / hbar_c for ich in range(l_Nch)])
        
        for jch in range(l_Nch):
            if (a_PotFunc_NLO is a_Params_NLO is None):
                l_TmpPhi = odeint(func_sch_ode,     l_InitPhi[:, jch], l_Range, 
                                  args=(a_PotFunc__LO, a_Params__LO, l_Evec, l_Mvec, a_l))[1]
            else:
                l_TmpPhi = odeint(func_sch_ode_NLO, l_InitPhi[:, jch], l_Range, 
                                  args=(a_PotFunc__LO, a_Params__LO, a_PotFunc_NLO, a_Params_NLO, 
                                        l_Evec, l_Mvec, a_l))[1]
            
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
