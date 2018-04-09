#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Sun May 28 00:27:38 JST 2017
### Brief : (1) plot Z-factor

import numpy as np

ifpathPS   = None
ifpathSS   = None
mass       = 1.0
Bsize      = 1

rcolumn    = 1

def main():
    CorrPS = read_correlators(ifpathPS)
    CorrSS = read_correlators(ifpathSS)
    
    CorrPS_jk = make_jk_samples(CorrPS[:,:,rcolumn], Bsize)
    if (CorrPS_jk is None):
        return
    CorrSS_jk = make_jk_samples(CorrSS[:,:,rcolumn], Bsize)
    if (CorrSS_jk is None):
        return
    
    plot_Zfac_jkerr(CorrPS_jk, CorrSS_jk, mass)

def read_correlators(a_path):
    from glob    import glob, iglob
    from os.path import dirname, basename
    
    rCorr = np.array([np.loadtxt(ifname) for ifname in iglob(a_path)]) #; print rCorr; quit()
    
    Nfile = len(rCorr[:,0,0])
    Ntime = len(rCorr[0,:,0])
    Ncolm = len(rCorr[0,0,:]) ; print("# Nfile=%d, Ntime=%d, Ncolm=%d" % (Nfile,Ntime,Ncolm))
    
    path_anti = dirname(a_path) + "/anti" + basename(a_path)
    
    if (True if (len(glob(path_anti)) == 0) else False):
        return rCorr
    else:
        print("# Read anti-hadron...")
        rCorr_anti = np.array([np.loadtxt(ifname) for ifname in iglob(path_anti)])
        rCorr_fb   = np.empty((Nfile, Ntime, Ncolm))
        for ifile in range(Nfile):
            for it in range(Ntime):
                for icol in range(Ncolm):
                    rCorr_fb[ifile,it,icol] = (rCorr     [ifile, it,icol] + 
                                               rCorr_anti[ifile,-it,icol]) / 2.0
        return rCorr_fb

def make_jk_samples(idata, a_Bsize):
    Ndata = len(idata[:,0])
    
    if (Ndata <= a_Bsize):
        print("\nERROR (Ndata=%d): Unexpected #.file "
              "(#.data <= #.bin), exit.\n" % Ndata); return None
    if (a_Bsize != 0):
        if (Ndata % a_Bsize != 0):
            print("\nERROR (Ndata=%d): Unexpected #.file "
                  "(#.data %% #.bin != 0), exit.\n" % Ndata); return None
    
    if (a_Bsize != 0):
        Nbin  = int(Ndata / a_Bsize)
        fdata = np.array([(np.sum(idata,axis=0)-np.sum(idata[i*a_Bsize:(i+1)*a_Bsize,:],axis=0))
                          for i in range(Nbin)]) / (Ndata-a_Bsize)
    else:
        fdata = idata
    
    return fdata

def plot_Zfac_jkerr(a_CorrPS, a_CorrSS, a_mass):
    Ndata = len(a_CorrPS[:,0])
    
    print("#")
    for it in range(len(a_CorrPS[0,:])):
        factor = np.sqrt(2.0*a_mass / np.exp(-a_mass*it))
        if ((a_CorrSS[:,it] <= 0).any()):
            ZfacP = 0
            ZfacS = 0
        else:
            ZfacP = a_CorrPS[:,it]/np.sqrt(a_CorrSS[:,it]) * factor
            ZfacS =                np.sqrt(a_CorrSS[:,it]) * factor 
        print("%d %1.16e %1.16e %1.16e %1.16e" % (it, 
                                                  np.mean(ZfacP), np.std(ZfacP) * np.sqrt(Ndata-1),
                                                  np.mean(ZfacS), np.std(ZfacS) * np.sqrt(Ndata-1)))
    print("#")

###### Functions for arguments
def usage(ARGV0):
    from os.path import basename
    
    print("\nusage: python %s '[ifpath.PS]' '[ifpath.SS]' {options}\n" % basename(ARGV0))
    print("option:")
    print("      -mass  [Mass (Latt.Unit)] Default = ", end=''); print(mass)
    print("      -bsize [bin-size        ] Default = ", end=''); print(Bsize)
    print

def set_args(ARGC, ARGV):
    global ifpathPS, ifpathSS, mass, Bsize
    from glob import glob
    
    if (ARGC < 3):
        usage(ARGV[0]); return -1
    if (ARGV[1][0] == '-' or ARGV[2][0] == '-'):
        usage(ARGV[0]); return -1
    
    ifpathPS = ARGV[1]
    ifpathSS = ARGV[2]
    NfilePS  = len(glob(ifpathPS))
    NfileSS  = len(glob(ifpathSS))
    if (NfilePS == 0):
        print("\nERROR: No such a file '%s', exit.\n" % ifpathPS); return -1
    if (NfileSS == 0):
        print("\nERROR: No such a file '%s', exit.\n" % ifpathSS); return -1
    if (NfilePS != NfileSS):
        print("\nERRPR: Different #.file in PS and SS, exit.\n"); return -1
    
    for i in range(3, ARGC):
        if (ARGV[i][0] == '-'):
            if   (ARGV[i] == '-mass'):
                mass = float(ARGV[i+1])
            elif (ARGV[i] == '-bsize'):
                Bsize = int(ARGV[i+1])
            else:
                print("\nERROR: Invalid argument '%s'" % ARGV[i])
                usage(ARGV[0]); return -1
    return 0

def check_args():
    print("# === Check Arguments ===")
    print("# if.PS = ", end=''); print(ifpathPS)
    print("# if.SS = ", end=''); print(ifpathSS)
    print("# mass  = ", end=''); print(mass)
    print("# Bsize = ", end=''); print(Bsize)
    print("# =======================")

if (__name__ == "__main__"):
    from sys import argv
    
    argc = len(argv)
    
    if (set_args(argc, argv) == 0):
        check_args() #; quit()
        main()
