# -*- coding: utf-8 -*-

"""The module to calculate the laplacian for a 3D-field on a lattice."""

from numpy import roll

def lap(a_data):
    """
    The function to calculate the laplacian for a 3D-field on a lattice.
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    return: r_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    Note: Assuming the periodic boundary condition.
    """
    r_data = -6.0 * a_data
    for ia in range(3):
        for pm1 in [-1, +1]:
            r_data += roll(a_data, pm1, axis=ia)
    return r_data

def lap_4th_proc(a_data):
    """
    The function to calculate the laplacian for a 3D-field on a lattice.
    ~~~ with 4th-order precision for 2nd-derivative ~~~
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    return: r_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    Note: Assuming the periodic boundary condition.
    """
    r_data = -90.0 * a_data
    for ia in range(3):
        for pm1 in [-1, +1]:
            r_data += 16.0 * roll(a_data, pm1, axis=ia)
        for pm2 in [-2, +2]:
            r_data -=        roll(a_data, pm2, axis=ia)
    return r_data / 12.0
