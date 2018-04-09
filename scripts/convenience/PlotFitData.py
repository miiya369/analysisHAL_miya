#!/usr/bin/python

if __name__ == "__main__":
    import sys
    import math
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 2):
        print "usage: %s [Data_fit]" % argv[0]
        quit()
    
    ifile = open(argv[1].strip(), 'rb')
    fsize = os.path.getsize(argv[1].strip())
    
    FitType = struct.unpack('<i', ifile.read(4))[0]
    Nconf   = struct.unpack('<i', ifile.read(4))[0]
    Ndata   = struct.unpack('<i', ifile.read(4))[0]
        
    print("# FitType=%d" % FitType)
    print("# N.conf =%d" % Nconf)
    print("# N.data =%d" % Ndata)
    
    if (FitType == 1):
        print("FitType = 1 has not been implemented.")
        quit(1)
    
    CheckFsize = 12 + 2*8*Ndata + 8*Nconf*Ndata
    if (fsize != CheckFsize):
        print "ERROR: This is may not be Miyamoto-Format data."
        quit()
    
    crd  = [0.0] *  Ndata
    data = [0.0] * (Ndata*Nconf)
    
    for n in range(Ndata):
        crd[n] = struct.unpack('<d', ifile.read(8))[0]
        tmp    = struct.unpack('<d', ifile.read(8))[0]
        for i in range(Nconf):
            data[n + Ndata*i] = struct.unpack('<d', ifile.read(8))[0]
    
    ifile.close()
    
    for i in range(Nconf):
        for n in range(Ndata):
            print("%lf %lf" % (crd[n], data[n + Ndata*i]))
        print("\n")
