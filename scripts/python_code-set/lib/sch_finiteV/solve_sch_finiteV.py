# -*- coding: utf-8 -*-

"""The module to solve the Schrodinger equation in a finite volume."""

from numpy               import empty
from scipy.sparse        import lil_matrix
from scipy.sparse.linalg import eigsh

def xyz(a_x, a_y, a_z, a_Lsize):
    return (a_x + a_Lsize * (a_y + a_Lsize * a_z))

def solve_sch_Fvol(a_V, a_mu, a_L, a_Nret, a_calc_evec = False):
    """
    The function to solve the Schrodinger equation in the finite volume.
    
    For arguments,
    - a_V[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    return: (array(eval) {, array(evec) if a_calc_evec == True})
    
    Note: a_V and a_mu should be in Lattice Unit.
    """
    
    l_Ham  = lil_matrix((a_L**3, a_L**3))
    nghbrs = empty(6)
    
    for x in range(a_L):
        for y in range(a_L):
            for z in range(a_L):
                ixyz = xyz(x,y,z,a_L)
                
                nghbrs[0] = xyz((x+1)    %a_L, y           , z           ,a_L)
                nghbrs[1] = xyz( x           ,(y+1)    %a_L, z           ,a_L)
                nghbrs[2] = xyz( x           , y           ,(z+1)    %a_L,a_L)
                nghbrs[3] = xyz((x+a_L-1)%a_L, y           , z           ,a_L)
                nghbrs[4] = xyz( x           ,(y+a_L-1)%a_L, z           ,a_L)
                nghbrs[5] = xyz( x           , y           ,(z+a_L-1)%a_L,a_L)
                
                l_Ham[ixyz,ixyz] = 6.0 / (2.0*a_mu) + a_V[z,y,x]
                for mu in range(6):
                    l_Ham[ixyz,nghbrs[mu]] = -1.0 / (2.0*a_mu)
    
    return eigsh(l_Ham, which='SA', k=a_Nret, return_eigenvectors=a_calc_evec)
