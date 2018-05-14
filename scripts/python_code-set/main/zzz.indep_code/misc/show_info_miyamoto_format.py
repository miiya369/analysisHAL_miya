#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    from struct  import unpack
    
    if (argc != 2):
        exit("usage: python %s [ifname]" % basename(argv[0]))
    
    ifname = argv[1].strip()
    
    with open(ifname, 'rb') as ifile:
        MgcNum = unpack('<i', ifile.read(4))[0]
        if   (MgcNum == 9190517):
            print("This is miyamoto-format file (ComplexField_BASE)")
            print("   - aSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - xSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - ySIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - zSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - tSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - bSIZE = %d" % unpack('<i', ifile.read(4))[0])
        elif (MgcNum == 94421393):
            print("This is miyamoto-format file (STATISTICS<ComplexField_BASE>)")
            print("   - aSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - xSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - ySIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - zSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - tSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - bSIZE = %d" % unpack('<i', ifile.read(4))[0])
            print("   - Ndata = %d" % unpack('<i', ifile.read(4))[0])
        elif (MgcNum == 19900518):
            print("This is miyamoto-format file (STATISTICS<ComplexField_XYZ> "+
                  "with 1/48 data size reduction)")
            print("   - #.conf        = %d" % unpack('<i', ifile.read(4))[0])
            print("   - #.data/conf   = %d" % unpack('<i', ifile.read(4))[0])
            print("   - bytes of data = %d" % unpack('<i', ifile.read(4))[0])
        else:
            exit("This file is not miyamoto-format file")
