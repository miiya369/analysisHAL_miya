# -*- coding: utf-8 -*-

"""The module to input/output T-matrix."""

from numpy   import array, empty
from struct  import pack, unpack
from os.path import getsize

MagicNumber_Tmat = 369018

def input_Tmatrix(a_iFname, verbose_flg = True):
    """
    The function to input the T-matrix.
    
    FORMAT: (1, int) Magic number (=369018)
    FORMAT: (2, int) #.conf
    FORMAT: (3, int) #.data
    FORMAT: (4, int) #.ch
    FORMAT: (5, float * #.data) value for x-coordination
    FORMAT: (6, float * 2 * #.conf * #.data * #.ch * #.ch) value for T-matrix
    FORMAT: All data are Little endian
    
    return: (x-coord[#.data], Tmat[#.conf, #.data, #.channel, #.channel])
    
    Note: Tmat is supposed to complex-type.
    """    
    
    with open(a_iFname, 'rb') as ifile:
        if (unpack('<i', ifile.read(4))[0] != MagicNumber_Tmat):
            print("\nERROR: This file is not Miyamoto-format T-matrix file, exit.\n")
            return (None, None)
        l_Nconf = unpack('<i', ifile.read(4))[0]
        l_Ndata = unpack('<i', ifile.read(4))[0]
        l_Nch   = unpack('<i', ifile.read(4))[0]
        
        Expected_Fsize = 4*4 + 8*l_Ndata + 16*l_Nconf*l_Ndata*l_Nch*l_Nch
        if (getsize(a_iFname) != Expected_Fsize):
            print("\nERROR: Unexpected file size, exit.\n")
            return (None, None)
        
        r_DataCrd = array([unpack('<d', ifile.read(8))[0] for idata in range(l_Ndata)])
        r_Tmatrix = empty((l_Nconf, l_Ndata, l_Nch, l_Nch), dtype=complex)
        for iconf in range(l_Nconf):
            for idata in range(l_Ndata):
                for ich in range(l_Nch):
                    for jch in range(l_Nch):
                        tmp_real = unpack('<d', ifile.read(8))[0]
                        tmp_imag = unpack('<d', ifile.read(8))[0]
                        r_Tmatrix[iconf,idata,ich,jch] = complex(tmp_real, tmp_imag)
    
    if (verbose_flg):
        print("# Successful to input T-matrix data")
        print("# N.conf =", len(r_Tmatrix[:,0,0,0]))
        print("# N.data =", len(r_Tmatrix[0,:,0,0]))
        print("# N.ch   =", len(r_Tmatrix[0,0,:,0]))
    
    return (r_DataCrd, r_Tmatrix)

def output_Tmatrix(a_oFname, a_xData, a_Tmat, verbose_flg = True):
    """
    The function to output the T-matrix.
    
    FORMAT: (1, int) Magic number (=369018)
    FORMAT: (2, int) #.conf
    FORMAT: (3, int) #.data
    FORMAT: (4, int) #.ch
    FORMAT: (5, float * #.data) value for x-coordination
    FORMAT: (6, float * 2 * #.conf * #.data * #.ch * #.ch) value for T-matrix
    FORMAT: All data are Little endian
    
    For arguments,
    - a_xData[#.data]                               (1-dim ndarray)
    - a_Tmat [#.conf, #.data, #.channel, #.channel] (4-dim ndarray)
    
    Note1: Tmat is supposed to complex-type.
    Note2: #.conf, #.data and #.channel are got from a_Tmat
    """
    
    l_Nconf = len(a_Tmat[:,0,0,0])
    l_Ndata = len(a_Tmat[0,:,0,0])
    l_Nch   = len(a_Tmat[0,0,:,0])
    
    if (len(a_xData) != l_Ndata):
        print("\nWARNING: Different #.data between T-matrix and xData, Stop output.\n")
        return
    
    with open(a_oFname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNumber_Tmat))
        ofile.write(pack('<i', l_Nconf))
        ofile.write(pack('<i', l_Ndata))
        ofile.write(pack('<i', l_Nch  ))
        ofile.write(pack('<%dd' % l_Ndata, *a_xData))
        for iconf in range(l_Nconf):
            for idata in range(l_Ndata):
                for ich in range(l_Nch):
                    for jch in range(l_Nch):
                        ofile.write(pack('<d', a_Tmat[iconf,idata,ich,jch].real))
                        ofile.write(pack('<d', a_Tmat[iconf,idata,ich,jch].imag))
    
    if (verbose_flg):
        print("# Successful to output T-matrix data")
