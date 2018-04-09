#!/usr/bin/python

if __name__ == "__main__":
    from sys     import argv; argc = len(argv)
    from os.path import basename
    import numpy as np
    
    if (argc < 3):
        print("usage: python %s [ofile] [ifile1] [ifile2] ..." % basename(argv[0])); quit()
    
    opath    = argv[1].strip() #; print(opath)
    ipath    = argv[2:]        #; print(ipath)
    Naverage = len(ipath)      #; print(Naverage); quit()
    
    corr_ave = np.loadtxt(ipath[0])
    Ncolm    = len(corr_ave[0,:]) #; print(Ncolm); quit()
    for i in range(1, Naverage):
        corr_ave += np.loadtxt(ipath[i])
    corr_ave /= float(Naverage)
    
    np.savetxt(opath, corr_ave, fmt='%4d'+' %1.16e'*(Ncolm-1))
