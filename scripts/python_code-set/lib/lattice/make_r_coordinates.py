# -*- coding: utf-8 -*-

"""The module to make a r-binning data."""

from numpy import roll

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
                    o_r.append(np.sqrt(x**2+y**2+z**2))
                    o_d.append(idata[z,y,x])
    else:
        l_Lz2 = len(a_data[:,0,0]) // 2
        l_Ly2 = len(a_data[0,:,0]) // 2
        l_Lx2 = len(a_data[0,0,:]) // 2
        for z in range(-l_Lz2, l_Lz2):
            for y in range(-l_Ly2, l_Ly2):
                for x in range(-l_Lx2, l_Lx2):
                    o_r.append(np.sqrt(x**2+y**2+z**2))
                    o_d.append(idata[z,y,x])
    return np.array((o_r, o_d))

def make_rbin(a_data, a_bsize = 1):
    """
    The function to make a r-binning data.
    
    For arguments,
    - a_data[2, #.data] (2-dim ndarray)
    ( a_data[0] = r-coordinates, a_data[1] = values )
    
    return: r_data[3, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = values, r_data[2] = errors )
    """
    if (len(idata[0,:]) != len(idata[1,:])):
        print("\nERROR: #.data of r-coordinates and values must be the same.\n")
        return None
    o_r = []
    o_d = []
    oer = []
    for r in np.arange(0.0, np.max(idata[0,:]), a_bsize):
        idcs = np.where((r <= idata[0,:]) & (idata[0,:] < r+a_bsize))[0]
        Nb   = len(idcs)
        if (Nb != 0):
            br  = np.sum(idata[0,idcs]   ) / float(Nb)
            bd  = np.sum(idata[1,idcs]   ) / float(Nb)
            bd2 = np.sum(idata[1,idcs]**2) / float(Nb)
            o_r.append(br)
            o_d.append(bd)
            if (Nb == 1):
                oer.append(0.0)
            else:
                oer.append(np.sqrt(bd2-bd**2) / np.sqrt(Nb-1))
    return np.array((o_r, o_d, oer))
