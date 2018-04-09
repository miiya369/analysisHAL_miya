#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename, getsize
    from struct  import pack, unpack
    
    if (argc != 4):
        exit("usage: python %s [binary file] " % basename(argv[0]) +
             "[endian convert (0=no, 1=yes)] [component number (0,1,2,...)]")
    
    ifname     = argv[1].strip()
    endian_cnv = int(argv[2])
    index      = int(argv[3])
    
    if (endian_cnv != 0 and endian_cnv != 1):
        exit("usage: python %s [binary file] " % basename(argv[0]) +
             "[endian convert (0=no, 1=yes)] [component number (0,1,2,...)]")
    if (index < 0):
        exit("usage: python %s [binary file] " % basename(argv[0]) +
             "[endian convert (0=no, 1=yes)] [component number (0,1,2,...)]")
    
    ifsize = getsize(ifname)
    if (ifsize % 8 != 0):
        exit("Unexpected file size: fsize % sizeof(double) != 0")
    dfsize = int(ifsize / 8)
    
    if (dfsize-1 < index):
        exit("Too large component: %d is max index" % (dfsize-1))
    
    with open(ifname, 'rb') as ifile:
        ifile.seek(8*index)
        tmp_d = unpack('d', ifile.read(8))[0]
        if (endian_cnv == 1):
            tmp_d = unpack(">d", pack("<d",tmp_d))[0]
        print("ifile[%d] = %1.16e" % (index, tmp_d))
