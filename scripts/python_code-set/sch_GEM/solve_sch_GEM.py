# -*- coding: utf-8 -*-

"""The module to solve the Schrodinger equation by using Gauss expansion method."""

from numpy        import sqrt, pi, array
from numpy.linalg import norm
from scipy.linalg import eigh

from sch_GEM.set_mat_GEM import *

def solve_sch_GEM(a_Nret, a_Range, a_Params, a_mu, a_fname, a_Np = 0, a_EigValOnly = True):
    """
    The function to solve the Schrodinger equation by using Gauss expansion method.
    
    For arguments,
    - a_Range [#.base]  (1-dim ndarray)
    - a_Params[#.param] (1-dim ndarray)
    
    return: (EigVal, EigVec)
    - EigVal[#.return]           (1-dim ndarray)
    - EigVec[#.base  , #.return] (2-dim ndarray)
    
    Note: If a_EigValOnly = True, return (EigVal, None)
    """
    l_Ham, l_Phi = set_mat(a_Range, a_Params, a_mu, a_fname, a_Np)
    l_Ret        = eigh(l_Ham, l_Phi, eigvals_only=a_EigValOnly, eigvals=(0, a_Nret-1))
    
    if (a_EigValOnly):
        return (l_Ret, None)
    else:
        ### For normalize
        #l_Ret[1][:,:] = array([[ l_Ret[1][ibase,iret] / norm(l_Ret[1][:,iret]) 
        #                         for iret in range(a_Nret)] 
        #                       for ibase in range(len(l_Ret[1][:,0]))])
        return (l_Ret[0], l_Ret[1])

def calc_from_HamMat(a_Nret, a_Range, a_HamMat, a_Np = 0, a_EigVal_ONLY = True):
    """
    The function to solve the Schrodinger equation by using Gauss expansion method.
    ~~~ From Hamiltonian matrixes ~~~
    
    For arguments,
    - a_Range[#.base]                 (1-dim array)
    - a_Ham  [#.base, #.base, #.conf] (3-dim array)
    
    return: (EigVal, EigVec)
    - EigVal[#.return, #.conf]           (2-dim array)
    - EigVec[#.base  , #.return, #.conf] (3-dim array)
    
    Note2: If a_EigVal_ONLY = True, returns are (EigVal, None)
    """
    
    l_Nbase = len(a_Range)
    l_Nconf = len(a_HamMat[0, 0, :])
    
    l_Phi = empty((l_Nbase, l_Nbase))
    for i in range(l_Nbase):
        for j in range(i, l_Nbase):
            nu_ij1      = a_Range[i]    * a_Range[j]
            nu_ij2      = a_Range[i]**2 + a_Range[j]**2
            l_Phi[i, j] = pow(2.0 * nu_ij1 / nu_ij2, 1.5)
            l_Phi[j, i] = l_Phi[i, j]
    
    r_val = empty((a_Nret, l_Nconf))
    if (a_EigVal_ONLY):
        r_vec = None
    else:
        r_vec = empty((l_Nbase, a_Nret, l_Nconf))
    
    print("#")
    for iconf in range(l_Nconf):
        if (isinstance(a_Np, int)): # For Coulomb potential
            l_Ham_C = add_Coulomb(a_Range, a_HamMat[:, :, iconf], a_Np)
        else:
            l_Ham_C = add_Coulomb_fpot_WS(a_Range, a_HamMat[:, :, iconf], a_Np[0], a_Np[1:])
        
        l_Ret = eigh(l_Ham_C, l_Phi, eigvals_only=a_EigVal_ONLY, eigvals=(0, a_Nret-1))
        
        if (a_EigVal_ONLY):
            for iret in range(a_Nret):
                r_val[iret, iconf] = l_Ret[iret]
        
        else:
            for iret in range(a_Nret):
                r_val[iret, iconf] = l_Ret[0][iret]
                for ibase in range(l_Nbase):
                    r_vec[ibase, iret, iconf] = l_Ret[1][ibase, iret]
        
        print("# conf = %d: End" % iconf)
    
    return (r_val, r_vec)
