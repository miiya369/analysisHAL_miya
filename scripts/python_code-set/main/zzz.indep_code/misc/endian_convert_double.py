#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename, getsize
    from struct  import pack, unpack
    
    if (argc != 3):
        exit("usage: python %s [ifname] [ofname]" % basename(argv[0]))
    
    ifname = argv[1].strip()
    ofname = argv[2].strip()
    
    ifsize = getsize(ifname)
    if (ifsize % 8 != 0):
        exit("Unexpected file size: fsize % sizeof(double) != 0")
    
    ifile = open(ifname, 'rb')
    ofile = open(ofname, 'wb')
    
    for i in range(int(ifsize / 8)):
        idata  =  unpack('d', ifile.read(8))[0]
        ofile.write(pack('d', unpack(">d", pack("<d", idata))[0]))

    ifile.close()
    ofile.close()
