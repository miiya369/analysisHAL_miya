#!/usr/bin/python

from __future__ import print_function

import numpy as np
from struct import pack, unpack

### =========================== Main =========================== ###

### Choose the data which will be made jackknife samples
def main (iconf_lst, oconf_lst, ifpath, ofpath):
    ifpaths = get_paths(iconf_lst, ifpath)
    ofpaths = get_paths(oconf_lst, ofpath)
    if (ifpaths is None or ofpaths is None):
        return -1
    Nconf = len(ifpaths) #; print(Nconf)
    Nbin  = len(ofpaths) #; print(Nbin); return 0
    if (Nconf % Nbin != 0):
        print("\nERROR: Unexpected #.file, #.conf(=%d) %% #.bin(=%d) != 0.\n" % (Nconf, Nbin)); return -1
    
    #return jkbin_textdata(ofpaths, ifpaths)
    #return jkbin_bindata (ofpaths, ifpaths)
    #return jkbin_wave48  (ofpaths, ifpaths)
    return -1

### ============================================================ ###

def get_paths (fconf_lst, fpath):
    if (not 'REPCONF' in fpath):
        print("\nERROR: Put 'REPCONF' in i/ofpath, which will be replaced by conf name in i/oconf.lst\n")
        return None
    with open(fconf_lst, 'r') as fconf:
        conf_lists = fconf.readlines() #; print(conf_lists)
        Nconf      = len(conf_lists)
    return [fpath.replace('REPCONF', conf_lists[i].strip()) for i in range(Nconf)]

def jkbin_textdata (ofpaths, ifpaths):
    idata   = np.array([np.loadtxt(ifpaths[i]) for i in range(len(ifpaths))])
    idata2  = np.reshape(idata,  (len(idata[:,0,0]), len(idata[0,:,0])* len(idata[0,0,:])))
    jkdata  = make_jk_samples(idata2, len(ofpaths))
    if (jkdata is None):
        return -1
    jkdata2 = np.reshape(jkdata, (len(jkdata[:,0]), len(idata[0,:,0]), len(idata[0,0,:])))
    for i in range(len(ofpaths)):
        np.savetxt(ofpaths[i], jkdata2[i,:,:], fmt='%4d'+' %1.16e'*(len(jkdata2[i,0,:])-1))
    return 0

def jkbin_bindata (ofpaths, ifpaths):
    idata  = np.array([np.fromfile(ifpaths[i], '>d') for i in range(len(ifpaths))])
    jkdata = make_jk_samples(idata, len(ofpaths))
    if (jkdata is None):
        return -1
    for i in range(len(ofpaths)):
        with open(ofpaths[i], 'wb') as ofdata:
            ofdata.write(pack('>%dd' % len(jkdata[i,:]), *jkdata[i,:]))
    return 0

def jkbin_wave48 (ofpaths, ifpaths):
    header, length = read_wave48_head(ifpaths[0]) #; print(length); quit()
    if (header is length is None):
        return -1
    idata  = np.array([read_wave48_body(ifpaths[i], length) for i in range(len(ifpaths))])
    if (np.any(idata == None)):
        return -1
    jkdata = make_jk_samples(idata, len(ofpaths))
    if (jkdata is None):
        return -1
    for i in range(len(ofpaths)):
        write_wave48(ofpaths[i], header, length, jkdata[i,:])
    return 0

def make_jk_samples (idata, Nbin):
    Nconf = len(idata[:,0])
    if (Nconf % Nbin != 0):
        print("\nERROR: Unexpected #.file, #.conf(=%d) %% #.bin(=%d) != 0.\n" % (Nconf, Nbin))
        return None
    Bsize = Nconf // Nbin
    return np.array([(np.sum(idata, axis=0)-np.sum(idata[i*Bsize:(i+1)*Bsize, :], axis=0))
                     for i in range(Nbin)]) / float(Nconf-Bsize)

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
    
    if (argc != 5):
        exit("usage: python %s [iconf.lst] [oconf.lst] [ifpath] [ofpath]" % basename(argv[0]) + "\n"
             "     :" + "\n"
             "     : Put 'REPCONF' in i/ofpath, which will be replaced by conf name in i/oconf.lst")
    
    iconf_lst = argv[1]
    oconf_lst = argv[2]
    ifpath    = argv[3]
    ofpath    = argv[4]
    
    if(main(iconf_lst, oconf_lst, ifpath, ofpath) != 0):
        exit("ERROR EXIT.")
