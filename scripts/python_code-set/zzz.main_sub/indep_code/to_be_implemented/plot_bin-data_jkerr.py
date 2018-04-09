#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Sun May 28 00:27:38 JST 2017
### Brief : (1) Input bin-data (file size should be 2*Lsize^3 * 8 bytes)
### Brief : (2) Estimate jack-knife error
### Brief : (3) Output for gnuplot

import sys
import numpy as np

N_additional_option = 1

def main(ARGC, ARGV):
    ### Input data & Error check
    Bsize = int(ARGV[1].strip())
    idata = np.array([np.fromfile(ARGV[i], '>d') for i in range(1 + N_additional_option, ARGC)])
    Ndata = len(idata[:,0])
    Fsize = len(idata[0,:])
    Lsize = int((Fsize/2)**(1.0/3.0)+0.001)
    if (not check_args(ARGC, Bsize, idata, Ndata, Fsize, Lsize)):
        return False
    
    ### Make jack-knife sample (Skip if Bsize == 0)
    if (Bsize != 0):
        Nbin  = Ndata / Bsize
        fdata = np.array([(np.sum(idata, axis=0)-np.sum(idata[i*Bsize:(i+1)*Bsize, :], axis=0)) 
                          for i in range(Nbin)]) / (Ndata-Bsize)
    else:
        Nbin  = Ndata
        fdata = idata
    
    ### Additional calculation for data
    #for i in range(Nbin):
    #    tmp_d = A1_projection(fdata[i, :], Lsize)
    #    for n in range(Fsize):
    #        fdata[i, n] = tmp_d[n]
    
    ### Output with jack-knife error
    factor = np.sqrt(Nbin-1)
    
    print("# R | Re(mean) | Re(error) | Im(mean) | Im(error)")
    for x in range(Lsize/2):
        for y in range(Lsize/2):
            for z in range(Lsize/2):
                #if (x > y > z):
                #    continue
                R      = np.sqrt(x**2+y**2+z**2)
                ixyz   = x + Lsize*(y + Lsize*z)
                mean_r = np.mean(fdata[:, 0 + 2*ixyz])
                mean_i = np.mean(fdata[:, 1 + 2*ixyz])
                err_r  = np.std (fdata[:, 0 + 2*ixyz]) * factor
                err_i  = np.std (fdata[:, 1 + 2*ixyz]) * factor
                print("%lf %e %e %e %e" % (R, mean_r, err_r, mean_i, err_i))
    
    return True

def check_args(a_argc, a_Bsize, a_data, a_Ndata, a_Fsize, a_Lsize):
    if (a_Ndata != a_argc-1 - N_additional_option):
        print("ERROR (Ndata=%d): Failed to read files, exit." % a_Ndata)
        return False
    if (a_Fsize != 2*a_Lsize**3):
        print("ERROR (Fsize=%d,Lsize=%d): Unexpected file size (not 2*Lsize^3 * 8 bytes), exit." % (a_Fsize,a_Lsize))
        return False
    for i in range(a_Ndata):
        if (len(a_data[i,:]) != a_Fsize):
            print("ERROR (Fsize=%d): Unexpected file size (different file size), exit." % len(a_data[i,:]))
            return False
    if (a_Bsize < 0):
        print("DEBUG MODE: Check Input")
        print("   - Ndata = %d" % a_Ndata)
        print("   - Fsize = %d" % a_Fsize)
        print("   - Lsize = %d" % a_Lsize)
        return False
    if (a_Ndata <= a_Bsize):
        print("ERROR (Ndata=%d): Unexpected #.file (#.data <= #.bin), exit." % a_Ndata)
        return False
    if (a_Bsize != 0):
        if (a_Ndata % a_Bsize != 0):
            print("ERROR (Ndata=%d): Unexpected #.file (#.data %% #.bin != 0), exit." % a_Ndata)
            return False
    return True

#def A1_projection(a_data, a_Lsize):

if (__name__ == "__main__"):
    argv = sys.argv; argc = len(argv)
    if (argc < 2 + N_additional_option):
        print("usage: python %s [bin-size] [bin-data(s)]..." % argv[0]); quit()
    
    main(argc, argv)
