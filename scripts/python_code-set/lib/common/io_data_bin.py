# -*- coding: utf-8 -*-

"""The module to input/output binary data."""

MagicNum_BinData_Miya = 19900518
MagicNum_BinData_NXYZ = 94421393

from numpy  import fromfile, array
from struct import pack

def input_bin_data(a_ifname, calc_err = True, verbose_flg = True):
    """
    The function to input the (miyamoto format) binary data.
    
    return: (yData[#.conf, #.data], xData[#.data] {, eData[#.data] for calc_err = True})
    """
    with open(a_ifname, 'rb') as ifile:
        l_MagicNum = fromfile(ifile, '<i4', 1)[0]
        
        if (l_MagicNum != MagicNum_BinData_Miya):
            print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum)
            return (None, None, None)
        
        l_Nconf = fromfile(ifile, '<i4', 1)[0]
        l_Ndata = fromfile(ifile, '<i4', 1)[0]
        l_NByte = fromfile(ifile, '<i4', 1)[0]
        
        if   (l_NByte == 4):
            r_xData = fromfile(ifile, '<i4' ,         l_Ndata)
            r_yData = fromfile(ifile, '<i4' , l_Nconf*l_Ndata).reshape((l_Nconf, l_Ndata))
        elif (l_NByte == 8):
            r_xData = fromfile(ifile, '<f8' ,         l_Ndata)
            r_yData = fromfile(ifile, '<f8' , l_Nconf*l_Ndata).reshape((l_Nconf, l_Ndata))
        elif (l_NByte == 16):
            r_xData = fromfile(ifile, '<f8' ,         l_Ndata)
            r_yData = fromfile(ifile, '<c16', l_Nconf*l_Ndata).reshape((l_Nconf, l_Ndata))
        else:
            print("\nERROR: Cannot input with #.byte =", l_NByte, ".\n")
            return None
    
    if (verbose_flg):
        print("# Successful to input Binary data")
        print("# N.conf = %d" % l_Nconf)
        print("# N.data = %d" % l_Ndata)
    
    if (calc_err):
        return (r_yData, r_xData, std(r_yData, axis=0) * sqrt(l_Nconf-1))
    else:
        return (r_yData, r_xData)

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
    
    if   (a_xData.dtype == 'int32'   and a_yData.dtype == 'int32'     ):
        l_NByte = 4
    elif (a_xData.dtype == 'float64' and a_yData.dtype == 'float64'   ):
        l_NByte = 8
    elif (a_xData.dtype == 'float64' and a_yData.dtype == 'complex128'):
        l_NByte = 16
    else:
        print("\nERROR: Cannot output with dtype (xData,yData)=(",
              a_xData.dtype,",",a_yData.dtype,").\n")
        return
    
    with open(a_ofname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNum_BinData_Miya))
        ofile.write(pack('<i', l_Nconf))
        ofile.write(pack('<i', l_Ndata))
        ofile.write(pack('<i', l_NByte))
        
        if   (l_NByte == 4):
            ofile.write(pack('<%di' %         l_Ndata  , *a_xData))
            ofile.write(pack('<%di' % l_Nconf*l_Ndata  , *a_yData.flatten()))
        elif (l_NByte == 8):
            ofile.write(pack('<%dd' %         l_Ndata  , *a_xData))
            ofile.write(pack('<%dd' % l_Nconf*l_Ndata  , *a_yData.flatten()))
        elif (l_NByte == 16):
            ofile.write(pack('<%dd' %         l_Ndata  , *a_xData))
            ofile.write(pack('<%dd' % l_Nconf*l_Ndata*2,
                             *array((a_yData.flatten().real,
                                     a_yData.flatten().imag)).T.flatten()))
    if (verbose_flg):
        print("# Successful to output bin data.")

def input_bin_data_nxyz(a_ifname, data_type = "complex", verbose_flg = True):
    """
    The function to input the (xyz-complex field w/all-conf) binary data.
    
    return: Data[#.conf, #.Zsite, #.Ysite, #.Xsite] (4-dim complex/float ndarray)
    """
    with open(a_ifname, 'rb') as ifile:
        l_MagicNum = fromfile(ifile, '<i4', 1)[0]
        
        if (l_MagicNum != MagicNum_BinData_NXYZ):
            print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum)
            return None
        
        l_Asize = fromfile(ifile, '<i4', 1)[0]
        l_Xsize = fromfile(ifile, '<i4', 1)[0]
        l_Ysize = fromfile(ifile, '<i4', 1)[0]
        l_Zsize = fromfile(ifile, '<i4', 1)[0]
        l_Tsize = fromfile(ifile, '<i4', 1)[0]
        l_Bsize = fromfile(ifile, '<i4', 1)[0]
        l_Nconf = fromfile(ifile, '<i4', 1)[0]
        
        if (l_Asize != 1 or l_Tsize != 1 or l_Bsize != 1):
            print("\nERROR: Input for (Asize,Tsize,Bsize) = (%d,%d,%d) is not allowed.\n" % 
                  (l_Asize, l_Tsize, l_Bsize))
            return None
        
        if   (data_type ==  "double"):
            r_Data = (fromfile(ifile, '<c16', l_Nconf*l_Zsize*l_Ysize*l_Xsize).real
                      .reshape((l_Nconf, l_Zsize, l_Ysize, l_Xsize)))
        elif (data_type == "complex"):
            r_Data = (fromfile(ifile, '<c16', l_Nconf*l_Zsize*l_Ysize*l_Xsize)
                      .reshape((l_Nconf, l_Zsize, l_Ysize, l_Xsize)))
        else:
            print("\nERROR: Only the data types 'complex' or 'double' can be specified.\n" +
                  "Input the (xyz-complex field w/all-conf) binary data was failed.\n")
            return None
        
    if (verbose_flg):
        print("# Successful to input Binary data")
        print("# X-size = %d" % l_Xsize)
        print("# Y-size = %d" % l_Ysize)
        print("# Z-size = %d" % l_Zsize)
        print("# N.conf = %d" % l_Nconf)
    
    return r_Data

def output_bin_data_nxyz(a_ofname, a_data, verbose_flg = True):
    """
    The function to output the (xyz-complex field w/all-conf) binary data.
    
    For arguments,
    - a_data[#.conf, #.Zsite, #.Ysite, #.Xsite] (4-dim complex/float ndarray)
    
    Note: #.conf, #.Zsite, #.Ysite and #.Xsite are got from a_data
    """    
    l_Nconf = len(a_data[:,0,0,0])
    l_Zsize = len(a_data[0,:,0,0])
    l_Ysize = len(a_data[0,0,:,0])
    l_Xsize = len(a_data[0,0,0,:])
    
    with open(a_ofname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNum_BinData_NXYZ))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', l_Xsize))
        ofile.write(pack('<i', l_Ysize))
        ofile.write(pack('<i', l_Zsize))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', l_Nconf))
        
        ofile.write(pack('<%dd' % l_Nconf*l_Zsize*l_Ysize*l_Xsize * 2,
                         *array((a_data.flatten().real, 
                                 a_data.flatten().imag)).T.flatten()))
    if (verbose_flg):
        print("# Successful to output bin data.")

def input_bin_data_nt(a_ifname, data_type = "complex", verbose_flg = True):
    """
    The function to input the (t-complex field w/all-conf) binary data.
    
    return: Data[#.conf, #.Tsite] (2-dim complex/float ndarray)
    """
    with open(a_ifname, 'rb') as ifile:
        l_MagicNum = fromfile(ifile, '<i4', 1)[0]
        
        if (l_MagicNum != MagicNum_BinData_NXYZ):
            print("\nERROR: Invalid magic number '%d', exit.\n" % l_MagicNum)
            return None
        
        l_Asize = fromfile(ifile, '<i4', 1)[0]
        l_Xsize = fromfile(ifile, '<i4', 1)[0]
        l_Ysize = fromfile(ifile, '<i4', 1)[0]
        l_Zsize = fromfile(ifile, '<i4', 1)[0]
        l_Tsize = fromfile(ifile, '<i4', 1)[0]
        l_Bsize = fromfile(ifile, '<i4', 1)[0]
        l_Nconf = fromfile(ifile, '<i4', 1)[0]
        
        if (l_Asize != 1 or l_Xsize != 1 or l_Ysize != 1 or l_Zsize != 1 or l_Bsize != 1):
            print("\nERROR: Input for (Asize,Xsize,Ysize,Zsize,Bsize) = " + 
                  "(%d,%d,%d,%d,%d) is not allowed.\n" %
                  (l_Asize,l_Xsize,l_Ysize,l_Zsize,l_Bsize))
            return None
        
        if   (data_type ==  "double"):
            r_Data = (fromfile(ifile, '<c16', l_Nconf*l_Tsize).real
                      .reshape((l_Nconf, l_Tsize)))
        elif (data_type == "complex"):
            r_Data = (fromfile(ifile, '<c16', l_Nconf*l_Tsize)
                      .reshape((l_Nconf, l_Tsize)))
        else:
            print("\nERROR: Only the data types 'complex' or 'double' can be specified.\n" +
                  "Input the (xyz-complex field w/all-conf) binary data was failed.\n")
            return None
        
    if (verbose_flg):
        print("# Successful to input Binary data")
        print("# T-size = %d" % l_Tsize)
        print("# N.conf = %d" % l_Nconf)
    
    return r_Data

def output_bin_data_nt(a_ofname, a_data, verbose_flg = True):
    """
    The function to output the (t-complex field w/all-conf) binary data.
    
    For arguments,
    - a_data[#.conf, #.Tsite] (2-dim complex/float ndarray)
    
    Note: #.conf and #.Tsite are got from a_data
    """    
    l_Nconf = len(a_data[:,0])
    l_Tsize = len(a_data[0,:])
    
    with open(a_ofname, 'wb') as ofile:
        ofile.write(pack('<i', MagicNum_BinData_NXYZ))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', l_Tsize))
        ofile.write(pack('<i', 1))
        ofile.write(pack('<i', l_Nconf))
        
        ofile.write(pack('<%dd' % l_Nconf*l_Tsize * 2,
                         *array((a_data.flatten().real, 
                                 a_data.flatten().imag)).T.flatten()))
    if (verbose_flg):
        print("# Successful to output bin data.")
