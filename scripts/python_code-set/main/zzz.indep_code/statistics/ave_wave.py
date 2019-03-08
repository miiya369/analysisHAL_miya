#!/usr/bin/python

from __future__ import print_function

import os
import numpy as np
import glob
from struct import pack, unpack

### =========================== Main =========================== ###

### Choose the data which will be averaged
def main (idir, odir, clst):
    #ave_data = lambda ofpath, ifpaths: average_textdata(ofpath, ifpaths)
    ave_data = lambda of, ifs: average_bindata (of, ifs)
    #ave_data = lambda ofpath, ifpaths: average_wave48  (ofpath, ifpaths)
    
    for conf_in in open(clst, 'r'):
        conf_t = conf_in.strip().replace('Aconfig','tconfig')
        conf_x = conf_in.strip().replace('Aconfig','xconfig')
        conf_y = conf_in.strip().replace('Aconfig','yconfig')
        conf_z = conf_in.strip().replace('Aconfig','zconfig')
        
        ifnames_t = glob.glob(idir+"/"+conf_t+"/*"+conf_t+"*")
        ifnames_x = glob.glob(idir+"/"+conf_x+"/*"+conf_x+"*")
        ifnames_y = glob.glob(idir+"/"+conf_y+"/*"+conf_y+"*")
        ifnames_z = glob.glob(idir+"/"+conf_z+"/*"+conf_z+"*")
        
        if (len(ifnames_t) != len(ifnames_x) or
            len(ifnames_x) != len(ifnames_y) or
            len(ifnames_y) != len(ifnames_z)):
            print("ERROR: #.data in (t,x,y,z)-configs are differ, exit")
            return -1
        if (len(ifnames_t) == 0):
            print("ERROR: No such (t,x,y,z)-configs files, exit")
            return -1
        
        print("#.time of NBS wave = %d" % len(ifnames_t))
        
        if not os.path.exists(odir+"/"+conf_in.strip()):
            os.makedirs(odir+"/"+conf_in.strip())
        
        for i in range(len(ifnames_t)):
            ifiles = [ifnames_t[i], ifnames_x[i], ifnames_y[i], ifnames_z[i]]
            
            ifbase_t = os.path.basename(ifnames_t[i]).replace("tconfig","Aconfig")
            ifbase_x = os.path.basename(ifnames_x[i]).replace("xconfig","Aconfig")
            ifbase_y = os.path.basename(ifnames_y[i]).replace("yconfig","Aconfig")
            ifbase_z = os.path.basename(ifnames_z[i]).replace("zconfig","Aconfig")
            
            if (ifbase_t != ifbase_x or
                ifbase_x != ifbase_y or
                ifbase_y != ifbase_z):
                print("ERROR: Base name for (t,x,y,z)-configs are differ, exit")
                return -1
            
            ofile  = odir+"/"+conf_in.strip()+"/"+ifbase_t
            ave_data(ofile, ifiles)
        
    return 0

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
    
    if (argc != 4):
        exit("usage: python %s [in dir] [out dir] [conf lst]" % basename(argv[0]))
    
    if (main(argv[1].strip(), argv[2].strip(), argv[3].strip()) != 0):
        exit("ERROR EXIT.")
