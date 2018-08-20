# -*- coding: utf-8 -*-

"""The module to make a r-binning data."""

from numpy import any, append, argsort, sqrt, array, arange, max, where, sum, roll, dot
from scipy.linalg import svd, diagsvd
from lattice.rotation_projection import rot_proj
from misc_QM.special_functions   import sph_harm_A1_xyz

def make_r_coord(a_Lsize, A1_reduction = True):
    """
    The function to make the r-coordinate on the periodic cubic lattice.
    
    For arguments,
    - a_Lsize (scalar integer)
    - A1_reduction = True : Output only (x <= y <= z <= L/2)
    
    return: r_data[2, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = The number of points in the same r )    
    """
    o_r = array([])
    o_c = array([])
    if (A1_reduction):
        for z in range(a_Lsize//2+1):
            for y in range(z+1):
                for x in range(y+1):
                    r = sqrt(x**2+y**2+z**2)
                    if (any(o_r == r)):
                        o_c[o_r == r] += 1
                    else:
                        o_r = append(o_r, r)
                        o_c = append(o_c, 1)
    else:
        for z in range(-a_Lsize//2, a_Lsize//2):
            for y in range(-a_Lsize//2, a_Lsize//2):
                for x in range(-a_Lsize//2, a_Lsize//2):
                    r = sqrt(x**2+y**2+z**2)
                    if (any(o_r == r)):
                        o_c[o_r == r] += 1
                    else:
                        o_r = append(o_r, r)
                        o_c = append(o_c, 1)
    
    idx_new = argsort(o_r)
    return array((array(o_r)[idx_new], array(o_c)[idx_new]))

def xyzdata_to_rdata(a_data, A1_reduction = True):
    """
    The function to convert (periodic) xyz-data to r-data.
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - A1_reduction = True : Output only (x <= y <= z <= L/2)
    
    return: r_data[5, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = values,
    ..r_data[2] = z-coordinates, r_data[3] = y-coordinates,
    ..r_data[4] = x-coordinates)
    
    Note1: #.Lsite is got from a_data.
    Note2: When A1_reduction = True, assuming #.Zsite == #.Ysite == #.Xsite.
    """
    o_r = []
    o_d = []
    o_x = []
    o_y = []
    o_z = []
    if (A1_reduction):
        if (len(a_data[:,0,0]) != len(a_data[0,:,0]) or
            len(a_data[:,0,0]) != len(a_data[0,0,:])):
            print("\nERROR: #.site (Nx,Ny,Nz) must be the same for A1-reduction.\n")
            return None
        l_L = len(a_data[:,0,0])
        for z in range(l_L//2+1):
            for y in range(z+1):
                for x in range(y+1):
                    o_x.append(x)
                    o_y.append(y)
                    o_z.append(z)
                    o_r.append(sqrt(x**2+y**2+z**2))
                    o_d.append(a_data[z,y,x])
    else:
        l_Lz2 = len(a_data[:,0,0]) // 2
        l_Ly2 = len(a_data[0,:,0]) // 2
        l_Lx2 = len(a_data[0,0,:]) // 2
        for z in range(-l_Lz2, l_Lz2):
            for y in range(-l_Ly2, l_Ly2):
                for x in range(-l_Lx2, l_Lx2):
                    o_x.append(x)
                    o_y.append(y)
                    o_z.append(z)
                    o_r.append(sqrt(x**2+y**2+z**2))
                    o_d.append(a_data[z,y,x])
    return array((o_r, o_d, o_z, o_y, o_x))

def make_rbin(a_data, a_bsize = 1.0):
    """
    The function to make a r-binning data.
    
    For arguments,
    - a_data[2, #.data] (2-dim ndarray)
    ( a_data[0] = r-coordinates, a_data[1] = values )
    
    return: r_data[3, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, r_data[1] = values, r_data[2] = errors )
    """
    if (len(a_data[0]) != len(a_data[1])):
        print("\nERROR: #.data of r-coordinates and values must be the same.\n")
        return None
    o_r = []
    o_d = []
    oer = []
    for r in arange(0.0, max(a_data[0]), a_bsize):
        idcs = where((r <= a_data[0]) & (a_data[0] < r+a_bsize))[0]
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

def divide_wave_L04(a_data, a_crit = 0.1, up_to_L6 = False):
    """
    The function to divide a wave function into the w-funcs with l = 0 and 4.
    (If up_to_L6 = True, the w-funcs with l = 0, 4 and 6 are calculated.)
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    return: r_data[3 or 4, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, 
    ..r_data[1] = w-func with l=0, r_data[2] = w-func with l=4,
    ..r_data[3] = w-func with l=6 if up_to_L6 = True)
    
    Note: #.Lsize is got from len(a_data[0,0,:])
    """
    l_rcrd  = make_r_coord(len(a_data[0,0,:]))
    l_rdata = xyzdata_to_rdata(rot_proj(a_data.real))
    l_idx_r = [where(l_rdata[0]==l_rcrd[0,i])[0] for i in range(len(l_rcrd[0]))]
    if (any(array([len(l_idx_r[i]) for i in range(len(l_idx_r))]) != l_rcrd[1])):
        print("\nERROR: #.point at the same r is not match.\n")
        return None
    
    l_Y00 = sph_harm_A1_xyz(0,0, 1,1,1)
    orcrd = []
    opsi0 = []
    opsi4 = []
    opsi6 = []
    if (up_to_L6):
        Np_min = 3
        Matfnc = lambda x,y,z,iNp: array([[l_Y00,
                                           sph_harm_A1_xyz(4,0,x[i],y[i],z[i]),
                                           sph_harm_A1_xyz(6,0,x[i],y[i],z[i])] 
                                          for i in range(iNp)])
    else:
        Np_min = 2
        Matfnc = lambda x,y,z,iNp: array([[l_Y00,
                                           sph_harm_A1_xyz(4,0,x[i],y[i],z[i])] 
                                          for i in range(iNp)])
    for i in range(len(l_idx_r)):
        Np = len(l_idx_r[i])
        if (Np < Np_min):
            continue
        Psi = l_rdata[1,l_idx_r[i]]
        P_z = l_rdata[2,l_idx_r[i]]
        P_y = l_rdata[3,l_idx_r[i]]
        P_x = l_rdata[4,l_idx_r[i]]
        Mat = Matfnc(P_x,P_y,P_z, Np)
        U, S, Vh = svd(Mat)
        if (any(abs(S[-1]/S[0])<a_crit)):
            continue
        owave = dot(Vh.T, dot(diagsvd(1.0/S,len(Vh),len(U)), dot(U.T, Psi)))
        opsi0.append(owave[0])
        opsi4.append(owave[1])
        if (up_to_L6):
            opsi6.append(owave[2])
        orcrd.append(l_rcrd[0][i])
    
    if (up_to_L6):
        return array((orcrd,opsi0,opsi4,opsi6))
    else:
        return array((orcrd,opsi0,opsi4))

def rbin_L04_separation(a_data, a_bsize = 1.0, a_crit = 0.1, up_to_L6 = False):
    """
    The function to separate wave functions with l = 0 and 4.
    (If up_to_L6 = True, the w-funcs with l = 0, 4 and 6 are separated.)
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    return: r_data[7 or 10, #.data] (2-dim ndarray)
    ( r_data[0] = r-coordinates, 
    ..r_data[1] = g0 (w-func with l=0), r_data[2] = d/dr g0, r_data[3] = d2/dr2 g0,
    ..r_data[4] = g4 (w-func with l=4), r_data[5] = d/dr g4, r_data[6] = d2/dr2 g4
    .{if up_to_L6 = True,
    ..r_data[7] = g6 (w-func with l=6), r_data[8] = d/dr g6, r_data[9] = d2/dr2 g6})
    
    Note: #.Lsize is got from len(a_data[0,0,:])
    """
    l_rdata = xyzdata_to_rdata(rot_proj(a_data.real))
    
    l_Y00 = sph_harm_A1_xyz(0,0, 1,1,1)
    orcrd = []
    opsi0 = []; odpsi0 = []; od2psi0 = []
    opsi4 = []; odpsi4 = []; od2psi4 = []
    opsi6 = []; odpsi6 = []; od2psi6 = []
    if (up_to_L6):
        Np_min = 9
        Matfnc = lambda dr,x,y,z,iNp: array([[l_Y00, l_Y00 * dr[i], l_Y00 * dr[i]**2,
                                              sph_harm_A1_xyz(4,0,x[i],y[i],z[i]),
                                              sph_harm_A1_xyz(4,0,x[i],y[i],z[i]) * dr[i],
                                              sph_harm_A1_xyz(4,0,x[i],y[i],z[i]) * dr[i]**2,
                                              sph_harm_A1_xyz(6,0,x[i],y[i],z[i]),
                                              sph_harm_A1_xyz(6,0,x[i],y[i],z[i]) * dr[i],
                                              sph_harm_A1_xyz(6,0,x[i],y[i],z[i]) * dr[i]**2]
                                             for i in range(iNp)])
    else:
        Np_min = 6
        Matfnc = lambda dr,x,y,z,iNp: array([[l_Y00, l_Y00 * dr[i], l_Y00 * dr[i]**2,
                                              sph_harm_A1_xyz(4,0,x[i],y[i],z[i]),
                                              sph_harm_A1_xyz(4,0,x[i],y[i],z[i]) * dr[i],
                                              sph_harm_A1_xyz(4,0,x[i],y[i],z[i]) * dr[i]**2]
                                             for i in range(iNp)])
    
    for r in arange(0.0, max(l_rdata[0]), a_bsize):
        r0   = r + a_bsize / 2.0
        idcs = where((r <= l_rdata[0]) & (l_rdata[0] < r+a_bsize))[0]
        Np   = len(idcs)
        if (Np < Np_min):
            continue
        dr  = l_rdata[0,idcs] - r0
        Psi = l_rdata[1,idcs]
        P_z = l_rdata[2,idcs]
        P_y = l_rdata[3,idcs]
        P_x = l_rdata[4,idcs]
        Mat = Matfnc(dr,P_x,P_y,P_z, Np)
        U, S, Vh = svd(Mat)
        if (any(abs(S[-1]/S[0])<a_crit)):
            continue
        owave = dot(Vh.T, dot(diagsvd(1.0/S,len(Vh),len(U)), dot(U.T, Psi)))
        opsi0.append(owave[0]); odpsi0.append(owave[1]); od2psi0.append(owave[2])
        opsi4.append(owave[3]); odpsi4.append(owave[4]); od2psi4.append(owave[5])
        if (up_to_L6):
            opsi6.append(owave[6]); odpsi6.append(owave[7]); od2psi6.append(owave[8])
        orcrd.append(r0)
    
    if (up_to_L6):
        return array((orcrd,opsi0,odpsi0,od2psi0,opsi4,odpsi4,od2psi4,opsi6,odpsi6,od2psi6))
    else:
        return array((orcrd,opsi0,odpsi0,od2psi0,opsi4,odpsi4,od2psi4))
