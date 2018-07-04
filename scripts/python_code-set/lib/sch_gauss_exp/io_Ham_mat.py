# -*- coding: utf-8 -*-

"""The module to input/output Hamiltonian Matrixes."""

MagicNum_HamMat = 91905178

from numpy  import array, empty
from struct import pack, unpack

def input_Ham_mat(a_ifname, verbose_flg = True):
    """
    The function to input the Hamiltonian Matrixes.
    
    return: (max_r, range[#.base], Ham[#.conf, #.base, #.base])
    """
    
    with open(a_ifname, 'rb') as ifile:
        l_MagicNum = unpack('<i', ifile.read(4))[0]
        
        if (l_MagicNum != MagicNum_HamMat):
            print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum)
            return (None, None, None, None)
        
        l_Nbase   = unpack('<i', ifile.read(4))[0]
        l_Nconf   = unpack('<i', ifile.read(4))[0]
        l_max_r   = unpack('<d', ifile.read(8))[0]
        
        r_range = array([unpack('<d', ifile.read(8))[0] for ibase in range(l_Nbase)])
        r_Ham   = empty((l_Nconf, l_Nbase, l_Nbase))
        
        for iconf in range(l_Nconf):
            for i in range(l_Nbase):
                for j in range(i, l_Nbase):
                    r_Ham[iconf,i,j] = unpack('<d', ifile.read(8))[0]
                    r_Ham[iconf,j,i] = r_Ham[iconf,i,j]
    
    if (verbose_flg):
        print("# Successful to input Hamiltonian matrixes")
        print("# N.base  =", l_Nbase)
        print("# N.conf  =", l_Nconf)
        print("# max r   =", l_max_r)
    
    return (l_max_r, r_range, r_Ham)

def output_Ham_mat(a_ofname, a_max_r, a_range, a_Ham, verbose_flg = True):
    """
    The function to output the Hamiltonian Matrixes.
    
    For arguments,
    - a_range[#.base]                 (1-dim ndarray)
    - a_Ham  [#.conf, #.base, #.base] (3-dim ndarray)
    
    Note: #.base and #.conf are got from a_Ham.
    """
    
    
    l_Nconf = len(a_Ham[:,0,0])
    l_Nbase = len(a_Ham[0,0,:])
    
    with open(a_ofname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNum_HamMat))
        ofile.write(pack('<i', l_Nbase))
        ofile.write(pack('<i', l_Nconf))
        ofile.write(pack('<d', a_max_r))
        
        ofile.write(pack('<%dd' % l_Nbase, *a_range))
        
        for iconf in range(l_Nconf):
            for i in range(l_Nbase):
                for j in range(i, l_Nbase):
                    ofile.write(pack('<d', a_Ham[iconf,i,j]))
    
    if (verbose_flg):
        print("# Successful to output Hamiltonian matrixes")
