#!/usr/bin/python

if __name__ == "__main__":
    import sys
    import math
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 3):
        print "usage: %s [Data_all] [ofile name]" % argv[0]
        quit()
    
    ifile = open(argv[1].strip(), 'rb')
    fsize = os.path.getsize(argv[1].strip())
    
    Nconf = struct.unpack('<i', ifile.read(4))[0]
    Ndata = struct.unpack('<i', ifile.read(4))[0]
    Bytes = struct.unpack('<i', ifile.read(4))[0]
    
    factor = Nconf - 1
    
    print("# N.conf=%d" % Nconf)
    print("# N.data=%d" % Ndata)
    
    CheckFsize = 12 + Bytes*Ndata + 16*Nconf*Ndata
    if (fsize != CheckFsize):
        print "ERROR: This is may not be Miyamoto-Format data."
        quit()
    
    crd  = [0.0] *  Ndata
    data = [0.0] * (Ndata*Nconf)
    
    for i in range(Ndata):
        if (Bytes == 4):
            crd[i] = struct.unpack('<i', ifile.read(4))[0]
        else:
            crd[i] = struct.unpack('<d', ifile.read(8))[0]
    
    for i in range(Nconf):
        for n in range(Ndata):
            data[n + Ndata*i] = struct.unpack('<d', ifile.read(8))[0]
            dummy             = struct.unpack('<d', ifile.read(8))[0]
    
    ifile.close()

    ofile   = open(argv[2].strip(), 'wb')
    FitType = 2
    
    ofile.write(struct.pack('<i', FitType))
    ofile.write(struct.pack('<i', Nconf))
    ofile.write(struct.pack('<i', Ndata))
    
    for n in range(Ndata):
        mean  = 0.0
        mean2 = 0.0
        for i in range(Nconf):
            mean  += data[n + Ndata*i]
            mean2 += data[n + Ndata*i] ** 2
        
        mean  /= float(Nconf);
        mean2 /= float(Nconf);
        err    = math.sqrt(factor * (mean2 - mean ** 2))
        
        ofile.write(struct.pack('<d', crd[n]))
        ofile.write(struct.pack('<d', err))
        for i in range(Nconf):
            ofile.write(struct.pack('<d', data[n + Ndata*i]))
            
    ofile.close()
