#!/usr/bin/python
#
# Author: Takaya Miyamoto
# Date  : Wed Nov  2 18:52:16 JST 2016
# Brief : Plot the compressed-NBS wave function
#

if __name__ == "__main__":
    import sys
    import math
    import cmath
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 3):
        print("usage: %s [NBSwave] [Spin= 0 or 1]" % argv[0])
        quit()
    
    fname = argv[1].strip()
    Spin  = argv[2].strip()
    
    if (Spin != "0" and Spin != "1"):
        print("usage: %s [NBSwave] [Spin= 0 or 1]" % argv[0])
        quit()
    
    spin = int(Spin)
    
### Read data
    ifile = open(fname, 'rb')
    
    dummy = struct.unpack('<i', ifile.read(4))[0]
    if (dummy != 16):
        print("ERROR: The file '%s' may be not Ishii-san's compressed NBS data." % fname)
        quit()
    
    Xsites = struct.unpack('<i', ifile.read(4))[0]
    Ysites = struct.unpack('<i', ifile.read(4))[0]
    Zsites = struct.unpack('<i', ifile.read(4))[0]
    length = struct.unpack('<i', ifile.read(4))[0]
    dummy  = ifile.read(4*2)[0]

    R = [0] * length
    for i in range(length):
        ixyz  = struct.unpack('<i', ifile.read(4))[0]
        ix    = ixyz % Xsites
        ixyz /= Xsites
        iy    = ixyz % Ysites
        ixyz /= Ysites
        iz    = ixyz % Zsites
        R[i]  = math.sqrt(ix*ix + iy*iy + iz*iz)
    
    dummy  = ifile.read(4*9)[0]

    tmpW = [0.0] * 16
    
    for i in range(length):
        for j in range(16):
            tmpr    = struct.unpack('<d', ifile.read(8))[0]
            tmpi    = struct.unpack('<d', ifile.read(8))[0]
            tmpW[j] = complex(tmpr, tmpi)
        
        if (spin == 0): # Spin projection (Not-normalize)
            Wave = tmpW[5] + tmpW[10] - tmpW[6] - tmpW[9]
        else:
            Wave = tmpW[5] + tmpW[10] + tmpW[6] + tmpW[9]
        
        print("%d %e %e" % (R[i], Wave.real, Wave.imag))
