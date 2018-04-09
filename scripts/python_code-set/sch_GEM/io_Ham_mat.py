# -*- coding: utf-8 -*-

"""The module to input/output Hamiltonian Matrixes."""

MagicNum_HamMat = 91905178

from numpy import empty
from numpy import append

def input_HamMat(a_iFname):
    """
    The function to input the Hamiltonian Matrixes.
    
    return: (range_a, max_r, Range[#.base], Ham[#.base, #.base, #.conf])
    """
    
    from struct import unpack
    
    ifile = open(a_iFname, 'rb')
    
    l_MagicNum = unpack('<i', ifile.read(4))[0]
    
    if (l_MagicNum != MagicNum_HamMat):
        print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum); ifile.close()
        return (None, None, None, None)
    
    l_Nbase   = unpack('<i', ifile.read(4))[0]
    l_Nconf   = unpack('<i', ifile.read(4))[0]
    l_range_a = unpack('<d', ifile.read(8))[0]
    l_max_r   = unpack('<d', ifile.read(8))[0]
    
    r_Range = empty((l_Nbase))
    r_Ham   = empty((l_Nbase, l_Nbase, l_Nconf))
    
    for i in range(l_Nbase):
        r_Range[i] = unpack('<d', ifile.read(8))[0]
    
    for iconf in range(l_Nconf):
        for i in range(l_Nbase):
            for j in range(i, l_Nbase):
                r_Ham[i, j, iconf] = unpack('<d', ifile.read(8))[0]
                r_Ham[j, i, iconf] = r_Ham[i, j, iconf]
    
    ifile.close()
    
    print("# Successful to input Hamiltonian matrixes")
    print("# N.base  ="),; print l_Nbase
    print("# N.conf  ="),; print l_Nconf
    print("# range a ="),; print l_range_a
    print("# max r   ="),; print l_max_r
    
    return (l_range_a, l_max_r, r_Range, r_Ham)

def output_HamMat(a_oFname, a_range_a, a_max_r, a_Range, a_Ham):
    """
    The function to output the Hamiltonian Matrixes.
    
    For arguments,
    - a_Range[#.base] (1-dim array)
    - a_Ham  [#.base, #.base, #.conf] (3-dim array)
    
    Note: #.base and #.conf are got from a_Ham.
    """
    
    from struct import pack
    
    l_Nbase = len(a_Ham[:, 0, 0])
    l_Nconf = len(a_Ham[0, 0, :])
    
    ofile = open(a_oFname, 'wb')
    
    ofile.write(pack('<i', MagicNum_HamMat))
    ofile.write(pack('<i', l_Nbase))
    ofile.write(pack('<i', l_Nconf))
    ofile.write(pack('<d', a_range_a))
    ofile.write(pack('<d', a_max_r))
    
    for i in range(l_Nbase):
        ofile.write(pack('<d', a_Range[i]))
    
    for iconf in range(l_Nconf):
        for i in range(l_Nbase):
            for j in range(i, l_Nbase):
                ofile.write(pack('<d', a_Ham[i, j, iconf]))
    
    ofile.close()
    print("# Successful to output Hamiltonian matrixes")
