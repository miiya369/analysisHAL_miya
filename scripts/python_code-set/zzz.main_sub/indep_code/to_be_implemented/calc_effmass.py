#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Sun May 28 00:27:38 JST 2017
### Brief : (1) plot effective mass
### Brief : (2) fit  effective mass

import numpy as np

ifpath     = None
t_min      = 6
t_max      = 12
Bsize      = 1

rcolumn    = 1
Init_param = (1.0, 0.01)

def main():
    Correlators = read_correlators(ifpath)
    
    Corr_jk = make_jk_samples(Correlators[:,:,rcolumn], Bsize)
    if (Corr_jk is None):
        return
    
    Nbin   = len(Corr_jk[:,0])
    Ntime  = len(Corr_jk[0,:])
    factor = np.sqrt(Nbin-1)
    
    plot_effmass_jkerr(Corr_jk)
    
    Corr_jk_err = np.array([np.std(Corr_jk[:,i])*factor for i in range(Ntime)])
    
    Effmass = np.array([fit_effmass(Corr_jk[i,:], Corr_jk_err, t_min, t_max) for i in range(Nbin)])
    
    print("# Effective mass = %15.16f +/- %15.16f\n#" % (np.mean(Effmass), np.std(Effmass)*factor))

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

def plot_effmass_jkerr(a_Corr):
    Ndata = len(a_Corr[:,0])
    
    print("#")
    for it in range(len(a_Corr[0,:])-1):
        if ((a_Corr[:,it] == 0).any()):
            Effmass = np.zeros(Ndata)
        else:
            if ((a_Corr[:,it+1]/a_Corr[:,it] <= 0).any()):
                Effmass = np.zeros(Ndata)
            else:
                Effmass = -np.log(a_Corr[:,it+1]/a_Corr[:,it])
        
        print("%d %1.16e %1.16e" % (it, np.mean(Effmass), np.std(Effmass) * np.sqrt(Ndata-1)))
    print("#")

def fit_effmass(a_Corr, a_Corr_err, a_tmin, a_tmax):
    from scipy.optimize import curve_fit
    
    return curve_fit(lambda it, Coe, Mas: Coe * np.exp(-Mas * it), 
                     np.array([i for i in range(a_tmin, a_tmax+1)]), 
                     a_Corr[a_tmin:a_tmax+1], sigma=a_Corr_err[a_tmin:a_tmax+1],
                     p0=np.array(Init_param))[0][1]

###### Functions for arguments
def usage(ARGV0):
    from os.path import basename
    
    print("\nusage: python %s '[ifpath]' {options}\n" % basename(ARGV0))
    print("option:")
    print("      -tmin  [Minimum fit range] Default = ", end=''); print(t_min)
    print("      -tmax  [Maximum fit range] Default = ", end=''); print(t_max)
    print("      -bsize [bin-size         ] Default = ", end=''); print(Bsize)
    print

def set_args(ARGC, ARGV):
    global ifpath, t_min, t_max, Bsize
    from glob import glob
    
    if (ARGC == 1):
        usage(ARGV[0]); return -1
    if (ARGV[1][0] == '-'):
        usage(ARGV[0]); return -1
    
    ifpath = ARGV[1]
    Nfile  = len(glob(ifpath))
    if (Nfile == 0):
        print("\nERROR: No such a file '%s', exit.\n" % ifpath); return -1
    
    for i in range(2, ARGC):
        if (ARGV[i][0] == '-'):
            if   (ARGV[i] == '-tmin'):
                t_min = int(ARGV[i+1])
            elif (ARGV[i] == '-tmax'):
                t_max = int(ARGV[i+1])
            elif (ARGV[i] == '-bsize'):
                Bsize = int(ARGV[i+1])
            else:
                print("\nERROR: Invalid argument '%s'" % ARGV[i])
                usage(ARGV[0]); return -1
    return 0

def check_args():
    print("# === Check Arguments ===")
    print("# ifile = ", end=''); print(ifpath)
    print("# t min = ", end=''); print(t_min)
    print("# t max = ", end=''); print(t_max)
    print("# Bsize = ", end=''); print(Bsize)
    print("# =======================")

if (__name__ == "__main__"):
    from sys import argv
    
    argc = len(argv)
    
    if (set_args(argc, argv) == 0):
        check_args() #; quit()
        main()
