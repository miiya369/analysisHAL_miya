# -*- coding: utf-8 -*-

"""The module to make a r-binning data."""

from numpy import any, append, sqrt, array, arange, max, where, sum

def make_r_coord(a_Lsize):
    """
    The function to make the r-coordinate on the periodic cubic lattice.
    
    For arguments,
    - a_Lsize (scalar integer)
    
    return: r_data[2, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = The number of points in the same r )    
    """
    o_r = array([])
    o_c = array([])
    for z in range(a_Lsize//2+1):
        for y in range(z+1):
            for x in range(y+1):
                r = sqrt(x**2+y**2+z**2)
                if (any(o_r == r)):
                    o_c[o_r == r] += 1
                else:
                    o_r = append(o_r, r)
                    o_c = append(o_c, 1)
    return array((o_r, o_c))

def xyzdata_to_rdata(a_data, A1_reduction = True):
    """
    The function to convert (periodic) xyz-data to r-data.
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - A1_reduction = True : Output only (x <= y <= z <= L/2)
    
    return: r_data[2, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = values )
    
    Note1: #.Lsite is got from a_data.
    Note2: When A1_reduction = True, assuming #.Zsite == #.Ysite == #.Xsite.
    """
    o_r = []
    o_d = []
    if (A1_reduction):
        if (len(a_data[:,0,0]) != len(a_data[0,:,0]) or
            len(a_data[:,0,0]) != len(a_data[0,0,:])):
            print("\nERROR: #.site (Nx,Ny,Nz) must be the same for A1-reduction.\n")
            return None
        l_L = len(a_data[:,0,0])
        for z in range(l_L//2+1):
            for y in range(z+1):
                for x in range(y+1):
                    o_r.append(sqrt(x**2+y**2+z**2))
                    o_d.append(a_data[z,y,x])
    else:
        l_Lz2 = len(a_data[:,0,0]) // 2
        l_Ly2 = len(a_data[0,:,0]) // 2
        l_Lx2 = len(a_data[0,0,:]) // 2
        for z in range(-l_Lz2, l_Lz2):
            for y in range(-l_Ly2, l_Ly2):
                for x in range(-l_Lx2, l_Lx2):
                    o_r.append(sqrt(x**2+y**2+z**2))
                    o_d.append(a_data[z,y,x])
    return array((o_r, o_d))

def make_rbin(a_data, a_bsize = 1):
    """
    The function to make a r-binning data.
    
    For arguments,
    - a_data[2, #.data] (2-dim ndarray)
    ( a_data[0] = r-coordinates, a_data[1] = values )
    
    return: r_data[3, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = values, r_data[2] = errors )
    """
    if (len(a_data[0,:]) != len(a_data[1,:])):
        print("\nERROR: #.data of r-coordinates and values must be the same.\n")
        return None
    o_r = []
    o_d = []
    oer = []
    for r in arange(0.0, max(a_data[0,:]), a_bsize):
        idcs = where((r <= a_data[0,:]) & (a_data[0,:] < r+a_bsize))[0]
        Nb   = len(idcs)
        if (Nb != 0):
            br  = sum(a_data[0,idcs]   ) / float(Nb)
            bd  = sum(a_data[1,idcs]   ) / float(Nb)
            bd2 = sum(a_data[1,idcs]**2) / float(Nb)
            o_r.append(br)
            o_d.append(bd)
            if (Nb == 1):
                oer.append(0.0)
            else:
                oer.append(sqrt(bd2-bd**2) / sqrt(Nb-1))
    return array((o_r, o_d, oer))
