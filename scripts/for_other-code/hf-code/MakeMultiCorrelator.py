#!/usr/bin/python

### Auther: Takaya Miyamoto
### Date  : Mon Oct 10 18:19:30 JST 2016
###
### Brief : Make the 2B-correlator from the NBS wave functions

if (__name__ == "__main__"):
    import sys
    import os.path
    import struct
    import math
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 4):
        print("usage: %s [NBS wave (Not 1/48 compressed)] [Min time] [Max time]" % argv[0])
        print("     : At the NBS files, \"REPTIME\" will be replaced by the given time (Triple digits).")
        quit()
    
    iFname = argv[1]
    Tmin   = int(argv[2])
    Tmax   = int(argv[3])
    
    Fsize  = os.path.getsize(iFname.replace('REPTIME', ("%03d" % Tmin)))
    if (Fsize % 16 != 0):
        sys.exit("Invalid file size: fsize % (sizeof(double) * 2) != 0")
    Vsize  = Fsize / 16
    
    for it in range(Tmax+1):
        if (it < Tmin):
            print("%d 0.0 0.0" % it)
        else:
            ifwave = open(iFname.replace('REPTIME', ("%03d" % it)), 'rb')
            tmp_r  = 0.0
            tmp_i  = 0.0
            for n in range(Vsize):
                tmp_r += struct.unpack('>d', ifwave.read(8))[0]
                tmp_i += struct.unpack('>d', ifwave.read(8))[0]
            print("%d %1.16e %1.16e" % (it, tmp_r, tmp_i))
