#!/usr/bin/python

import numpy as np
from struct import pack, unpack
from sys    import exit

def main (argc, argv):
    opath    = argv[1].strip() #; print(opath)
    ipath    = argv[2:]        #; print(ipath)
    Naverage = len(ipath)      #; print(Naverage)
    
    data_hed, length = read_wave_head(ipath[0]) #; print(length); quit()
    data_ave         = read_wave_comp(ipath[0], length)
    
    for i in range(1, Naverage):
        data_tmp  = read_wave_comp(ipath[i], length)
        data_ave += data_tmp
    data_ave /= float(Naverage)
    
    write_wave_comp(opath, data_hed, length, data_ave)

def read_wave_head (ifname):
    with open(ifname, 'rb') as fdata:
        fdata.seek(16)
        length = np.fromfile(fdata, '<i', count = 1)[0] #; print(length); quit()
        fdata.seek(0)
        rhead  = np.fromfile(fdata, '<i', count = (16+length)) #; print(rhead); quit()
        if (rhead[-1] != 256*length):
            exit("The file '%s' may be not Ishii-san's compressed NBS data." % ifname)
    
    return (rhead, length)

def read_wave_comp (ifname, length):
    with open(ifname, 'rb') as fdata:
        fdata.seek((16+length) * 4)
        rwave = np.fromfile(fdata, '<d', count = 32*length)
        dummy = unpack('<i', fdata.read(4))[0]
        if (dummy != 256*length):
            exit("The file '%s' may be not Ishii-san's compressed NBS data." % ifname)
    
    return rwave

def write_wave_comp (ofname, headers, length, odata):
    with open(ofname, 'wb') as fdata:
        for n in range(len(headers)):
            fdata.write(pack('<i', headers[n]))
        for n in range(len(odata)):
            fdata.write(pack('<d', odata[n]))
        
        fdata.write(pack('<i', 256*length))

if __name__ == "__main__":
    from sys     import argv; argc = len(argv)
    from os.path import basename
    
    if (argc < 3):
        print("usage: python %s [ofile] [ifile1] [ifile2] ..." % basename(argv[0])); quit()
    
    main(argc, argv)
