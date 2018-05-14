# -*- coding: utf-8 -*-
# Author: Takaya Miyamoto
# E-mail: takaya.miyamoto@yukawa.kyoto-u.ac.jp
# Since : Sun Feb 18 21:31:38 JST 2018
# Brief1: The definitions of the function for input the MNIST data
# Brief2: The definitions of the output (print) image data

from __future__ import print_function

MagicNum_MNIST_Image = 2051
MagicNum_MNIST_Label = 2049

def input_mnist_data_image(a_ifname, a_N = None, verbose_flg = True):
    """
    The function to input the MNIST image data.
    
    @ Return
    - r_data[#.data, #.row, #.column] (3-dim ndarray, dtype = ubyte)
    """
    from numpy  import array, reshape
    from struct import unpack
    
    with open(a_ifname, 'rb') as ifile:
        byte_order = '<'
        l_MagicNum = unpack(byte_order+'i', ifile.read(4))[0]
        if (l_MagicNum != MagicNum_MNIST_Image):
            ifile.seek(0,0)
            byte_order = '>'
            l_MagicNum = unpack(byte_order+'i', ifile.read(4))[0]
            if (l_MagicNum != MagicNum_MNIST_Image):
                print("\nERROR: Invalid magic number '%d'.\n" % l_MagicNum)
                return None
        
        l_Ndat = unpack(byte_order+'i', ifile.read(4))[0]
        l_Nrow = unpack(byte_order+'i', ifile.read(4))[0]
        l_Ncol = unpack(byte_order+'i', ifile.read(4))[0]
        if (a_N is not None):
            l_Ndat = min(l_Ndat, a_N)
        
        r_data = array([[[unpack(byte_order+'B', ifile.read(1))[0]
                          for icol in range(l_Ncol)]
                         for irow in range(l_Nrow)]
                        for idat in range(l_Ndat)])
    
    if (verbose_flg):
        print("# Successful to input the MNIST image data")
        print("# N.data   = %d" % l_Ndat)
        print("# N.row    = %d" % l_Nrow)
        print("# N.column = %d" % l_Ncol)
    
    return r_data

def input_mnist_data_label(a_ifname, a_N = None, verbose_flg = True):
    """
    The function to input the MNIST label data.
    
    @ Return
    - r_data[#.data] (1-dim ndarray, dtype = ubyte)
    """
    from numpy  import array, reshape
    from struct import unpack
    
    with open(a_ifname, 'rb') as ifile:
        byte_order = '<'
        l_MagicNum = unpack(byte_order+'i', ifile.read(4))[0]
        if (l_MagicNum != MagicNum_MNIST_Label):
            ifile.seek(0,0)
            byte_order = '>'
            l_MagicNum = unpack(byte_order+'i', ifile.read(4))[0]
            if (l_MagicNum != MagicNum_MNIST_Label):
                print("\nERROR: Invalid magic number '%d'.\n" % l_MagicNum)
                return None
        
        l_Ndat = unpack(byte_order+'i', ifile.read(4))[0]
        if (a_N is not None):
            l_Ndat = min(l_Ndat, a_N)
        
        r_data = array([unpack(byte_order+'B', ifile.read(1))[0] for idat in range(l_Ndat)])
    
    if (verbose_flg):
        print("# Successful to input the MNIST label data")
        print("# N.data   = %d" % l_Ndat)
    
    return r_data

def print_image(a_idata, max_val = 256):
    """
    The function to print the image data.
    
    @ Arguments
    - a_idata[#.row, #.column] (2-dim ndarray, dtype = ubyte)
    """
    plot_char = ("  ", "..", "::", "++", "**", "OO", "##", "@@")
    
    if (max_val % len(plot_char) != 0):
        print("\nERROR: max_val %% len(plot_char) != 0, cannot plot the image.\n")
    else:
        N = max_val // len(plot_char)
        for     irow in range(len(a_idata[:,0])):
            for icol in range(len(a_idata[0,:])):
                print(plot_char[a_idata[irow,icol]//N], end="")
            print("")
