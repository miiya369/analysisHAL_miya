#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from struct     import pack, unpack
    from subprocess import call
    
    otmp1 = 1
    
    with open("endian_check_tmp.bin", 'wb') as ofile:
        ofile.write(pack('i', otmp1))
    with open("endian_check_tmp.bin", 'rb') as ifile:
        itmp1 = unpack('<i', ifile.read(4))[0]
    
    call("rm endian_check_tmp.bin".split())
    
    if (otmp1 == itmp1):
        print("This machine is little endian.")
    else:
        print("This machine is big endian.")
