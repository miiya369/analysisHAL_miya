#!/usr/bin/python

if __name__ == "__main__":
    import sys
    import math
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 2):
        print "usage: %s [Data_all]" % argv[0]
        quit()
    
    ifile = open(argv[1].strip(), 'rb')
    fsize = os.path.getsize(argv[1].strip())
    
    Nconf = struct.unpack('<i', ifile.read(4))[0]
    Ndata = struct.unpack('<i', ifile.read(4))[0]
    Bytes = struct.unpack('<i', ifile.read(4))[0]
    
    factor = Nconf - 1
    
    print("# N.conf=%d" % Nconf)
    print("# N.data=%d" % Ndata)
    print("# coordinate | mean (real) | error (real) | mean (imag) | error (imag)")
    
    CheckFsize = 12 + Bytes*Ndata + 16*Nconf*Ndata
    if (fsize != CheckFsize):
        print "ERROR: This is may not be Miyamoto-Format data."
        quit()
    
    crd  = [0] *  Ndata
    data = [0] * (Ndata*Nconf*2)
    
    for i in range(Ndata):
        if (Bytes == 4):
            crd[i] = struct.unpack('<i', ifile.read(4))[0]
        else:
            crd[i] = struct.unpack('<d', ifile.read(8))[0]
    
    for i in range(Nconf):
        for j in range(Ndata):
            for k in range(2):
                data[k + 2*(j + Ndata*i)] = struct.unpack('<d', ifile.read(8))[0]
    
    ifile.close()
    
    for j in range(Ndata):
        meanR  = 0.0
        meanR2 = 0.0
        meanI  = 0.0
        meanI2 = 0.0
        for i in range(Nconf):
            meanR  += data[0 + 2*(j + Ndata*i)]
            meanR2 += data[0 + 2*(j + Ndata*i)] ** 2
            meanI  += data[1 + 2*(j + Ndata*i)]
            meanI2 += data[1 + 2*(j + Ndata*i)] ** 2
        
        meanR  /= float(Nconf);
        meanR2 /= float(Nconf);
        meanI  /= float(Nconf);
        meanI2 /= float(Nconf);
        errR    = math.sqrt(factor * (meanR2 - meanR ** 2))
        errI    = math.sqrt(factor * (meanI2 - meanI ** 2))
        
        print("%f %e %e %e %e" % (crd[j], meanR, errR, meanI, errI))
