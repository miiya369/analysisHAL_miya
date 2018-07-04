# -*- coding: utf-8 -*-

"""The module to input/output binary data."""

MagicNum_BinData_Miya = 19900518
MagicNum_BinData_NXYZ = 94421393

from common.statistics import make_mean_err
from numpy             import array, reshape
from struct            import pack, unpack

def input_bin_data(a_ifname, verbose_flg = True):
    """
    The function to input the (miyamoto format) binary data.
    
    return: (yData[#.conf, #.data], xData[#.data], eData[#.data])
    """
    
    with open(a_ifname, 'rb') as ifile:
        l_MagicNum = unpack('<i', ifile.read(4))[0]
        
        if (l_MagicNum != MagicNum_BinData_Miya):
            print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum)
            return (None, None, None)
        
        l_Nconf = unpack('<i', ifile.read(4))[0]
        l_Ndata = unpack('<i', ifile.read(4))[0]
        l_NByte = unpack('<i', ifile.read(4))[0]
        
        if (l_NByte != 8):
            print("\nERROR: Bytes of data '%d' != 8: Have not implemented yet, exit.\n" % l_NByte)
            return (None, None, None)
        
        r_xData = array([ unpack('<d', ifile.read(8))[0] for idata in range(l_Ndata)])
        r_yData = array([[unpack('<d', ifile.read(8))[0] for idata in range(l_Ndata)]
                         for iconf in range(l_Nconf)])
    
    if (verbose_flg):
        print("# Successful to input Binary data")
        print("# N.conf = %d" % l_Nconf)
        print("# N.data = %d" % l_Ndata)
    
    r_eData = array([make_mean_err(r_yData[:,idata])[1] for idata in range(l_Ndata)])
    
    return (r_yData, r_xData, r_eData)

def output_bin_data(a_ofname, a_yData, a_xData, verbose_flg = True):
    """
    The function to output the binary data.
    
    For arguments,
    - a_yData[#.conf, #.data] (2-dim ndarray)
    - a_xData[#.data]         (1-dim ndarray)
    
    Note: #.conf and #.data are got from a_yData.
    """
    
    l_Nconf = len(a_yData[:,0])
    l_Ndata = len(a_yData[0,:])
    
    with open(a_ofname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNum_BinData_Miya))
        ofile.write(pack('<i', l_Nconf))
        ofile.write(pack('<i', l_Ndata))
        ofile.write(pack('<i', 8))
        
        ofile.write(pack('<%dd' % l_Ndata, *a_xData))
        
        for iconf in range(l_Nconf):
            ofile.write(pack('<%dd' % l_Ndata, *a_yData[iconf,:]))
    
    if (verbose_flg):
        print("# Successful to output bin data.")

def input_bin_data_nxyz(a_ifname, data_type = "complex", verbose_flg = True):
    """
    The function to input the (xyz-complex field w/all-conf) binary data.
    
    return: Data[#.conf, #.Zsite, #.Ysite, #.Xsite] (4-dim complex/float ndarray)
    """
    
    with open(a_ifname, 'rb') as ifile:
        l_MagicNum = unpack('<i', ifile.read(4))[0]
        
        if (l_MagicNum != MagicNum_BinData_NXYZ):
            print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum)
            return None
        
        l_Asize = unpack('<i', ifile.read(4))[0]
        l_Xsize = unpack('<i', ifile.read(4))[0]
        l_Ysize = unpack('<i', ifile.read(4))[0]
        l_Zsize = unpack('<i', ifile.read(4))[0]
        l_Tsize = unpack('<i', ifile.read(4))[0]
        l_Bsize = unpack('<i', ifile.read(4))[0]
        l_Ndata = unpack('<i', ifile.read(4))[0]
        
        if (l_Asize != 1 or l_Tsize != 1 or l_Bsize != 1):
            print("\nERROR: Input for (Asize,Tsize,Bsize) = (%d,%d,%d) is not allowed, exit.\n" % 
                  (l_Asize, l_Tsize, l_Bsize))
            return None
        
        l_Data = [unpack('<d', ifile.read(8))[0] for i in range(l_Xsize*l_Ysize*l_Zsize*l_Ndata*2)]
    
    if   (data_type ==   "float"):
        r_Data = reshape(array([l_Data[0+2*i] for i in range(l_Xsize*l_Ysize*l_Zsize*l_Ndata)]), 
                         (l_Ndata, l_Zsize, l_Ysize, l_Xsize))
    elif (data_type == "complex"):
        r_Data = reshape(array([complex(l_Data[0+2*i], l_Data[1+2*i]) for i in range(l_Xsize*l_Ysize*l_Zsize*l_Ndata)]), 
                         (l_Ndata, l_Zsize, l_Ysize, l_Xsize))
    else:
        print("\nERROR: Only the data types 'complex' or 'float' can be input.\n" +
              "Input the (xyz-complex field w/all-conf) binary data was failed.\n")
        return None
    
    if (verbose_flg):
        print("# Successful to input Binary data")
        print("# X-size = %d" % l_Xsize)
        print("# Y-size = %d" % l_Ysize)
        print("# Z-size = %d" % l_Zsize)
        print("# N.conf = %d" % l_Ndata)
    
    return r_Data

def output_bin_data_nxyz(a_ofname, a_Data, verbose_flg = True):
    """
    The function to output the (xyz-complex field w/all-conf) binary data.
    
    For arguments,
    - a_Data[#.conf, #.Zsite, #.Ysite, #.Xsite] (4-dim complex/float ndarray)
    
    Note: #.conf, #.Zsite, #.Ysite and #.Xsite are got from a_Data
    """
    
    l_Nconf = len(a_Data[:,0,0,0])
    l_Zsize = len(a_Data[0,:,0,0])
    l_Ysize = len(a_Data[0,0,:,0])
    l_Xsize = len(a_Data[0,0,0,:])
    
    with open(a_ofname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNum_BinData_NXYZ))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', l_Xsize))
        ofile.write(pack('<i', l_Ysize))
        ofile.write(pack('<i', l_Zsize))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', l_Nconf))
        
        for i in range(l_Nconf):
            for z in range(l_Zsize):
                for y in range(l_Zsize):
                    for x in range(l_Zsize):
                        if   (isinstance(a_Data[i,x,y,z],   float)):
                            ofile.write(pack('<d', a_Data[i,x,y,z]))
                            ofile.write(pack('<d', 0.0))
                        elif (isinstance(a_Data[i,x,y,z], complex)):
                            ofile.write(pack('<d', a_Data[i,x,y,z].real))
                            ofile.write(pack('<d', a_Data[i,x,y,z].imag))
                        else:
                            print("\nERROR: Only the data types 'complex' or 'float' can be output.\n" +
                                  "Output the (xyz-complex field w/all-conf) binary data was failed.\n")
                            return None
    
    if (verbose_flg):
        print("# Successful to output bin data.")
