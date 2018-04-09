#!/usr/bin/python

if __name__ == "__main__":
    import sys
    import math
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 5):
        print "usage: %s [ifile: Vc 1S0] [ifile: Vc 3S1] [ofile: V0] [ofile: Vsigma]" % argv[0]
        quit()
    
    ifile1 = open(argv[1].strip(), 'rb')
    ifile2 = open(argv[2].strip(), 'rb')
    fsize1 = os.path.getsize(argv[1].strip())
    fsize2 = os.path.getsize(argv[2].strip())
    
    if (fsize1 != fsize2):
        print "ERROR: Two files are different size."
        quit()
    
    Nconf = struct.unpack('<i', ifile1.read(4))[0]
    Ndata = struct.unpack('<i', ifile1.read(4))[0]
    Bytes = struct.unpack('<i', ifile1.read(4))[0]
    
    CheckFsize = 12 + Bytes*Ndata + 16*Nconf*Ndata
    if (fsize1 != CheckFsize):
        print "ERROR: This is may not be Miyamoto-Format data."
        quit()
    if (Bytes != 8):
        print "ERROR: This is may not be potential data."
        quit()
    if (struct.unpack('<i', ifile2.read(4))[0] != Nconf):
        print "ERROR: Two files are different Nconf."
        quit()
    if (struct.unpack('<i', ifile2.read(4))[0] != Ndata):
        print "ERROR: Two files are different Ndata."
        quit()
    if (struct.unpack('<i', ifile2.read(4))[0] != Bytes):
        print "ERROR: Two files are different Bytes."
        quit()
    
    print("#.conf = %d" % Nconf)
    print("#.data = %d" % Ndata)
    
    crd    = [0] *  Ndata
    idata1 = [0] * (Ndata*Nconf*2)
    idata2 = [0] * (Ndata*Nconf*2)
    
    for i in range(Ndata):
        crd[i] = struct.unpack('<d', ifile1.read(8))[0]
        if (     struct.unpack('<d', ifile2.read(8))[0] != crd[i]):
            print "ERROR: Two files are different coordinate."
            quit()
    
    for i in range(Nconf):
        for j in range(Ndata):
            for k in range(2):
                idata1[k + 2*(j + Ndata*i)] = struct.unpack('<d', ifile1.read(8))[0]
                idata2[k + 2*(j + Ndata*i)] = struct.unpack('<d', ifile2.read(8))[0]
    
    ifile1.close()
    ifile2.close()
    
    print("END reading, START calculation")
    
    odata1 = [0] * (Ndata*Nconf*2)
    odata2 = [0] * (Ndata*Nconf*2)
    
    for i in range(Nconf):
        for j in range(Ndata):
            for k in range(2):
                odata1[k+2*(j+Ndata*i)] = (idata1[k+2*(j+Ndata*i)] + 3.0 * idata2[k+2*(j+Ndata*i)]) / 4.0
                odata2[k+2*(j+Ndata*i)] = (idata2[k+2*(j+Ndata*i)] -       idata1[k+2*(j+Ndata*i)]) / 4.0
    
    print("END calculation, START writing")
    
    ofile1 = open(argv[3].strip(), 'wb')
    ofile2 = open(argv[4].strip(), 'wb')
    
    ofile1.write(struct.pack('<i', Nconf))
    ofile2.write(struct.pack('<i', Nconf))
    ofile1.write(struct.pack('<i', Ndata))
    ofile2.write(struct.pack('<i', Ndata))
    ofile1.write(struct.pack('<i', Bytes))
    ofile2.write(struct.pack('<i', Bytes))
    
    for i in range(Ndata):
        ofile1.write(struct.pack('<d', crd[i]))
        ofile2.write(struct.pack('<d', crd[i]))
    
    for i in range(Nconf):
        for j in range(Ndata):
            for k in range(2):
                ofile1.write(struct.pack('<d', odata1[k + 2*(j + Ndata*i)]))
                ofile2.write(struct.pack('<d', odata2[k + 2*(j + Ndata*i)]))
    
    ofile1.close()
    ofile2.close()

    print("END writing")
