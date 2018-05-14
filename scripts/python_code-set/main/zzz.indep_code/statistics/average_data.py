#!/usr/bin/python

from __future__ import print_function

import numpy as np
from struct import pack, unpack

### =========================== Main =========================== ###

### Choose the data which will be averaged
def main (ofpath, ifpaths):
    #return average_textdata(ofpath, ifpaths)
    #return average_bindata (ofpath, ifpaths)
    #return average_wave48  (ofpath, ifpaths)
    return -1

### ============================================================ ###

def average_textdata (ofpath, ifpaths):
    ave_data = np.mean(np.array([np.loadtxt(ifpaths[i]) for i in range(len(ifpaths))]), axis=0)
    np.savetxt(ofpath, ave_data, fmt='%4d'+' %1.16e'*(len(ave_data[0, :])-1))
    return 0

def average_bindata (ofpath, ifpaths):
    ave_data = np.mean(np.array([np.fromfile(ifpaths[i], '>d') for i in range(len(ifpaths))]), axis=0)
    with open(ofpath, 'wb') as ofdata:
        ofdata.write(pack('>%dd' % len(ave_data), *ave_data)) 
    return 0

def average_wave48 (ofpath, ifpaths):
    header, length = read_wave48_head(ifpaths[0]) #; print(length); quit()
    if (header is length is None):
        return -1
    write_wave48(ofpath, header, length, 
                 np.mean(np.array([read_wave48_body(ifpaths[i], length) for i in range(len(ifpaths))]), axis=0))
    return 0

def read_wave48_head (ifname):
    with open(ifname, 'rb') as fdata:
        fdata.seek(16); length = np.fromfile(fdata, '<i', count =           1)[0] #; print(length); quit()
        fdata.seek( 0); header = np.fromfile(fdata, '<i', count = (16+length))    #; print(header); quit()
        if (header[-1] != 256*length):
            print("\nREAD ERROR: The file '%s' may be not Ishii-san's compressed NBS data.\n" % ifname)
            return (None, None)
    return (header, length)

def read_wave48_body (ifname, length):
    with open(ifname, 'rb') as fdata:
        fdata.seek((16+length)*4); body = np.fromfile(fdata, '<d', count = 32*length)
        dummy = unpack('<i', fdata.read(4))[0]
        if (dummy != 256*length):
            print("\nREAD ERROR: The file '%s' may be not Ishii-san's compressed NBS data.\n" % ifname)
            return None
    return body

def write_wave48 (ofname, header, length, body):
    with open(ofname, 'wb') as fdata:
        fdata.write(pack('<%di' % len(header), *header))
        fdata.write(pack('<%dd' % len(body  ), *body  ))
        fdata.write(pack('<i', 256*length))

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc < 3):
        exit("usage: python %s [ofile] [ifile1] [ifile2] ..." % basename(argv[0]))
    
    ofpath  = argv[1].strip() #; print(ofpath)
    ifpaths = argv[2:]        #; print(ifpaths); print(len(ifpaths)); exit()
    
    if (main(ofpath, ifpaths) != 0):
        exit("ERROR EXIT.")
