#!/usr/bin/python

if __name__ == "__main__":
    from sys     import argv; argc = len(argv)
    from os.path import basename
    from struct  import pack, unpack
    import numpy as np
    
    if (argc < 3):
        print("usage: python %s [ofile] [ifile1] [ifile2] ..." % basename(argv[0])); quit()
    
    opath    = argv[1].strip() #; print(opath)
    ipath    = argv[2:]        #; print(ipath)
    Naverage = len(ipath)      #; print(Naverage); quit()
    
    data_ave = np.fromfile(ipath[0], '>d')
    for i in range(1, Naverage):
        data_ave += np.fromfile(ipath[i], '>d')
    data_ave /= float(Naverage)
    
    with open(opath, 'wb') as fdata:
        for n in range(len(data_ave)):
            fdata.write(pack('>d', data_ave[n]))
