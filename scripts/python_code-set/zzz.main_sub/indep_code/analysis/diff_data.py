#!/usr/bin/python

from __future__ import print_function

import numpy as np
from struct import unpack

### =========================== Main =========================== ###

def main (ifpath_1, ifpath_2, criterion):
    #idata_1, idata_2 = read_textdata(ifpath_1, ifpath_2)
    #idata_1, idata_2 = read_bindata (ifpath_1, ifpath_2)
    #idata_1, idata_2 = read_wave48  (ifpath_1, ifpath_2)
    idata_1 = idata_2 = None
    if (idata_1 is idata_2 is None):
        return -1
    
    if (np.all(idata_1 == idata_2)):
        print("Same."); return 0
    else:
        Ndata = len(idata_1)
        fdiff = np.empty(Ndata)
        for i in range(Ndata):
            if (idata_1[i] == 0.0 and idata_1[i] == idata_2[i]):
                fdiff[i] = 0.0
            else:
                fdiff[i] = abs(2.0*(idata_1[i]-idata_2[i])/(idata_1[i]+idata_2[i]))
        
        Ndiff = len(fdiff[fdiff > criterion])
        print("Diff = %1.16f +/- %1.16f" % (np.mean(fdiff), np.std(fdiff, ddof=1)))
        print("#.(diff > %1.16f) = %d (in the elements %d)" % (criterion, Ndiff, len(idata_2)))
        print("Maximam diff = %1.16f" % np.max(fdiff))
        return 0

### ============================================================ ###

def read_textdata (ifpath_1, ifpath_2):
    idata_1 = np.reshape(np.loadtxt(ifpath_1), -1)
    idata_2 = np.reshape(np.loadtxt(ifpath_2), -1)
    if (len(idata_1) != len(idata_2)):
        print("File sizes are different."); return (None, None)
    return (idata_1, idata_2)

def read_bindata (ifpath_1, ifpath_2):
    idata_1 = np.fromfile(ifpath_1, '>d')
    idata_2 = np.fromfile(ifpath_2, '>d')
    if (len(idata_1) != len(idata_2)):
        print("File sizes are different."); return (None, None)
    return (idata_1, idata_2)

def read_wave48 (ifpath_1, ifpath_2):
    data_hed_1, length_1 = read_wave_head(ifpath_1)
    if (data_hed_1 is length_1 is None):
        return (None, None)
    data_hed_2, length_2 = read_wave_head(ifpath_2)
    if (data_hed_2 is length_2 is None):
        return (None, None)
    if (np.any(data_hed_1 != data_hed_2)):
        print("Headers are different."); return (None, None)
    idata_1 = read_wave_comp(ifpath_1, length_1)
    if (idata_1 is None):
        return (None, None)
    idata_2 = read_wave_comp(ifpath_2, length_2)
    if (idata_2 is None):
        return (None, None)
    return (idata_1, idata_2)

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

### ============================================================ ###
### ============================================================ ###

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 4):
        exit("usage: python %s [ifile1] [ifile2] [criterion]" % basename(argv[0]))
    
    ifpath_1  = argv[1]
    ifpath_2  = argv[2]
    criterion = float(argv[3].strip())
    
    if(main(ifpath_1, ifpath_2, criterion) != 0):
        exit("ERROR EXIT.")
