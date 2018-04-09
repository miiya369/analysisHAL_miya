# -*- coding: utf-8 -*-

"""The module to solve the Schrodinger equation by using Gauss expansion method."""

from numpy.linalg import norm
from scipy.linalg import eigh
from math         import sqrt, pi

import re

from SchGaussExp.SetMatrix import *

def solve_sch_GEM(a_Nret, a_Range, a_Params, a_mass, FitType, a_Np = 0, a_EigVal_ONLY = True):
    """
    The function to solve the Schrodinger equation by using Gauss expansion method.
    
    For arguments,
    - a_Range [#.base]  (1-dim array)
    - a_Params[#.param] (1-dim array)
    
    return: (EigVal, EigVec)
    - EigVal[#.return]           (1-dim array)
    - EigVec[#.base  , #.return] (2-dim array)
    
    Note2: If a_EigVal_ONLY = True, returns are (EigVal, None)
    """
    
    if (re.match('^SW$', FitType) is not None):
        l_Ham, l_Phi = set_mat_idogata(a_Range, a_Params, a_mass)
        
        # For Coulomb potential
        l_Ham_C = add_Coulomb(a_Range, l_Ham, a_Np)
    
    else:
        if (isinstance(a_mass, float)):
            if   (re.match('^[1-9]G$', FitType) is not None):
                l_Ham, l_Phi = set_mat_gauss(a_Range, a_Params, a_mass)
                #l_Ham, l_Phi = set_mat(a_Range, a_Params, FitType, a_mass)
            
            else:
                l_Ham, l_Phi = set_mat(a_Range, a_Params, FitType, a_mass)
            
            # For Coulomb potential
            l_Ham_C = add_Coulomb(a_Range, l_Ham, a_Np)
            
        else:
            if   (re.match('^[1-9]G$', FitType) is not None):
                l_Ham, l_Phi = set_mat_gauss_fpot_WS(a_Range, a_Params, a_mass[0], a_mass[1:])
                #l_Ham, l_Phi = set_mat_fpot_WS(a_Range, a_Params, FitType, a_mass[0], a_mass[1:])
            
            else:
                l_Ham, l_Phi = set_mat_fpot_WS(a_Range, a_Params, FitType, a_mass[0], a_mass[1:])
            
            # For Coulomb potential of Folding potential with WS-type density
            l_Ham_C = add_Coulomb_fpot_WS(a_Range, l_Ham, a_Np, a_mass[1:])
    
    l_Ret = eigh(l_Ham_C, l_Phi, eigvals_only=a_EigVal_ONLY, eigvals=(0, a_Nret-1))
    
    if (a_EigVal_ONLY):
        return (l_Ret, None)
    
    else:
        #for iret in range(a_Nret): # For normalize
        #    tmp_norm = norm(l_Ret[1][:, iret])
        #    for ibase in range(len(l_Ret[1][:, iret])):
        #        l_Ret[1][ibase, iret] /= tmp_norm
        
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
