#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename, getsize
    from struct  import pack, unpack
    
    if (argc < 3):
        exit("usage: python %s [Binary file] [Format1] [Format2] ...\n" % basename(argv[0]) +
             "     : Format should be 'int' or 'double'\n" +
             "     : To skip the byte, put [bytes (Normal Num.)]")
    
    ifname = argv[1].strip()
    fsize  = getsize(ifname)
    print("Start: file size = %d bytes" % fsize)
    
    formats = argv[2:]
    Bcount  = 0
    with open(ifname, 'rb') as ifile:    
        for i in range(len(formats)):
            if   (formats[i].strip() == "int"):
                if (Bcount + 4 > fsize):
                    break
                print("from %10d bytes: ( int32) %d"     % (Bcount, unpack('i', ifile.read(4))[0]))
                Bcount += 4
            elif  (formats[i].strip() == "int_"):
                if (Bcount + 4 > fsize):
                    break
                print("from %10d bytes: ( int32) %d"     % (Bcount, unpack(">i", pack("<i", unpack('i', ifile.read(4))[0]))[0]))
                Bcount += 4
            elif (formats[i].strip() == "double"):
                if (Bcount + 8 > fsize):
                    break
                print("from %10d bytes: (double) %1.16e" % (Bcount, unpack('d', ifile.read(8))[0]))
                Bcount += 8
            elif (formats[i].strip() == "double_"):
                if (Bcount + 8 > fsize):
                    break
                print("from %10d bytes: (double) %1.16e" % (Bcount, unpack(">d", pack("<d", unpack('d', ifile.read(8))[0]))[0]))
                Bcount += 8
            elif (formats[i].strip().isdigit()):
                if (Bcount + int(argv[i+2]) > fsize):
                    break
                print("from %10d bytes:          SKIP %d bytes" % (Bcount, int(formats[i])))
                ifile.read(int(formats[i]))
                Bcount += int(formats[i])
            else:
                exit("Invalid format '%s'" % formats[i])
    
    print("Finish")
