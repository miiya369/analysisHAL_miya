# -*- coding: utf-8 -*-

"""The module to calculate a projection for angular momentum on a lattice."""

from numpy                   import array, zeros, dot
from lattice.rotation_matrix import irep, rot_mat, rot_chr

def rot_proj(a_Data, a_state = "A1"):
    """
    The function to calculate a projection for angular momentum on a lattice.
    
    For arguments,
    - a_Data[#.site, #.site, #.site] (3-dim ndarray)
    - a_state (= 'A1', 'A2', 'E', 'T1', 'T2')
    
    return: r_Data[#.site, #.site, #.site] (3-dim ndarray)
    
    Note: #.site is got from a_Data
    """
    
    if (len(a_Data[:,0,0]) != len(a_Data[0,:,0]) or
        len(a_Data[:,0,0]) != len(a_Data[0,0,:])):
        print("\nERROR: #.site (Nx,Ny,Nz) should be the same for rot_proj.\n")
        return None
    
    l_L = len(a_Data[:,0,0])
    
    if   (a_state == "A1" or a_state == "A2"):
        factor = 1.0 / 24.0
    elif (a_state ==  "E"):
        factor = 2.0 / 24.0
    elif (a_state == "T1" or a_state == "T2"):
        factor = 3.0 / 24.0
    else:
        print("\nERROR: Invalid representation type (%s)\n" % a_state)
        return None
    
    r_Data = zeros((l_L, l_L, l_L))
    
    for i in range(24):
        Rmat = rot_mat(l_L    , *irep(i)   )
        Rchr = rot_chr(a_state,  irep(i)[0])
        for z in range(l_L):
            for y in range(l_L):
                for x in range(l_L):
                    xyz_in = array([x,y,z,1])
                    xyz = dot(Rmat, xyz_in) % l_L
                    r_Data[xyz[2],xyz[1],xyz[0]] += Rchr * a_Data[z,y,x]
    
    return r_Data * factor
