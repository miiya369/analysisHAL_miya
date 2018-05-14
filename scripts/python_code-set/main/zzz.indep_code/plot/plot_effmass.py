#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc, patches

xlim_all = (0, 15)
#xlim_all = None

ylim_all = (1, 1.5)
#ylim_all = None

BinSize     = 1
read_column = 3

tmin = 6
tmax = 9

#lat_spacing = 0.1215
lat_spacing = None

fit_type = "exp"
#fit_type = "cosh"

label_location = 'upper left'

ofname = None
#ofname = "/Users/miiya/Desktop/effmass_kaon.pdf"

### =========================== Main =========================== ###

def main (ifpaths):
    Corr = read_correlators(ifpaths[0])
    for i in range(1, len(ifpaths)):
        Corr = np.append(Corr, read_correlators(ifpaths[i]), axis=0)
    Nfile = len(Corr[:,0,0])
    Ntime = len(Corr[0,:,0])
    Ncolm = len(Corr[0,0,:])
    print("# Read correlators: Nfile=%d, Ntime=%d, Ncolm=%d" % (Nfile, Ntime, Ncolm))
    
    if   (BinSize == 0):
        Nbin    = Nfile; 
        Corr_jk = Corr[:,:,read_column]
    elif (Nfile % BinSize != 0 or BinSize == Nfile):
        print("\nERROR: Unexpected binsize.\n"); return -1
    else:
        Nbin    = Nfile // BinSize
        Corr_jk = make_jk_samples(Corr[:,:,read_column], Nbin)
        if (Corr_jk is None):
            return -1
        print("# Make jackknife samples: Nbin=%d" % Nbin)
    
    if (lat_spacing is None):
        factor_MeV = 1.0
        unit_str   = "[Lattice Unit]"
    else:
        factor_MeV = 197.327053 / lat_spacing
        unit_str   = "[MeV]"
    
    if (fit_type == "exp"):
        Effmass = calc_effmass_exp(Corr_jk) * factor_MeV
        fitmass =  fit_effmass_exp(Corr_jk, tmin, tmax)
    elif (fit_type == "cosh"):
        Effmass = calc_effmass_csh(Corr_jk) * factor_MeV
        fitmass =  fit_effmass_csh(Corr_jk, tmin, tmax)
    else:
        print("\nERROR: Unknown fit type, %s.\n" % fit_type); return -1
    
    #print_data_err(Effmass)
    
    fitmass_ave = np.mean(fitmass[:,1] * factor_MeV)
    fitmass_err = np.std (fitmass[:,1] * factor_MeV) * np.sqrt(Nbin-1)
    print("# Effective mass = %15.16f +/- %15.16f %s" % (fitmass_ave, fitmass_err, unit_str))
    
    # --- plot --- #
    rc('text' , usetex=True)
    rc('font' ,**{'family'    : 'Times New Roman',
                  'weight'    : 'bold',
                  'size'      : 11})
    rc('xtick',**{'labelsize' : 12})
    rc('ytick',**{'labelsize' : 12})
    rc('axes' ,**{'labelsize' : 12})
    rc('axes' ,**{'linewidth' : 2})
    
    plt.ion(); fig = plt.figure()
    fig.patch.set_facecolor('white')
    axe = fig.add_axes((0.15, 0.15, 0.8, 0.8))
    axe.patch.set_facecolor('white')
    
    axe.errorbar(np.arange(len(Effmass[0,:])), np.mean(Effmass, axis=0), 
                 yerr=(np.std(Effmass, axis=0) * np.sqrt(Nbin-1)), fmt='none', #fmt=None,
                 marker=None, label='Effective mass data', lw=3, ecolor='red', capthick=3, capsize=3)
    
    axe.add_patch(patches.Rectangle((tmin, fitmass_ave-fitmass_err), tmax-tmin, 2*fitmass_err, 
                                     color='blue', alpha=0.5))
    axe.plot((tmin, tmax), (fitmass_ave, fitmass_ave), lw=3, color='blue',
             label='Fit value (t = %d-%d): %6f +/- %6f %s' % 
             (tmin, tmax, fitmass_ave, fitmass_err, unit_str))
    
    if (ylim_all is not None):
        axe.set_ylim(ylim_all)
    if (xlim_all is not None):
        axe.set_xlim(xlim_all)
    axe.legend(numpoints=1, loc=label_location) #, bbox_to_anchor=(0.0, 0.0))
    axe.set_ylabel('Effective mass %s' % unit_str)
    axe.set_xlabel('time slice ($t-t_0$)')
    axe.yaxis.set_label_coords(-0.12,0.5)
    axe.xaxis.set_label_coords(0.5,-0.1)
    #axe.axhline(0.0, lw=2, color='black')
    axe.grid(which='major',color='gray',linestyle='--')
    
    if (ofname is not None):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    else:
        plt.show(); print("# push Enter to end"); input('# ') #raw_input('# ') # For v.2
    
    return 0

### ============================================================ ###

def read_correlators(ifpath, auto_fb = True):
    from glob    import glob, iglob
    from os.path import dirname, basename
    Corr        = np.array([np.loadtxt(ifname) for ifname in iglob(ifpath)])
    ifpath_anti = dirname(ifpath) + "/anti" + basename(ifpath)
    if (len(glob(ifpath_anti)) == 0 or not auto_fb):
        return Corr
    else:
        print("# Read anti-hadron...")
        Corr_anti = np.array([np.loadtxt(ifname) for ifname in iglob(ifpath_anti)])
        return (Corr + np.roll(Corr_anti[:,::-1,:], 1, axis=1)) / 2.0

def make_jk_samples (idata, Nbin):
    Nconf = len(idata[:,0])
    if (Nconf % Nbin != 0):
        print("\nERROR: Unexpected #.file, #.conf(=%d) %% #.bin(=%d) != 0.\n" % (Nconf, Nbin))
        return None
    Bsize = Nconf // Nbin
    return np.array([(np.sum(idata, axis=0)-np.sum(idata[i*Bsize:(i+1)*Bsize, :], axis=0))
                     for i in range(Nbin)]) / float(Nconf-Bsize)

def calc_effmass_exp (iCorr):
    Nconf       = len(iCorr[:,0])
    Ntime       = len(iCorr[0,:])
    ret_effmass = np.empty((Nconf, Ntime-1))
    for it in range(Ntime-1):
        if (np.all(iCorr[:,it] != 0)):
            tmp_d = iCorr[:,it+1] / iCorr[:,it]
            if (np.all(tmp_d > 0)):
                ret_effmass[:,it] = -np.log(tmp_d)
            else:
                ret_effmass[:,it] = np.zeros(Nconf)
        else:
            ret_effmass[:,it] = np.zeros(Nconf)
    return ret_effmass

def calc_effmass_csh (iCorr):
    Nconf       = len(iCorr[:,0])
    Ntime       = len(iCorr[0,:])
    ret_effmass = np.empty((Nconf, Ntime-1))
    for it in range(Ntime-1):
        if (np.all(iCorr[:,it] != 0)):
            tmp_d = (iCorr[:,it+1] + iCorr[:,it-1]) / (2.*iCorr[:,it])
            if (np.all(tmp_d > 1)):
                ret_effmass[:,it] = np.arccosh(tmp_d)
            else:
                ret_effmass[:,it] = np.zeros(Nconf)
        else:
            ret_effmass[:,it] = np.zeros(Nconf)
    return ret_effmass

def print_data_err (idata):
    idata_ave = np.mean(idata, axis=0)
    idata_err = np.std (idata, axis=0) * np.sqrt(len(idata[:,0])-1)
    for i in range(len(idata_ave)):
        print("%d %15.16e %15.16e" % (i, idata_ave[i], idata_err[i]))

def fit_effmass_exp (iCorr, itmin, itmax):
    from scipy.optimize import curve_fit
    Nconf     = len(iCorr[:,0])
    iCorr_err = np.std(iCorr, axis=0) * np.sqrt(Nconf-1)
    return np.array([curve_fit(lambda it, Coe, Mas: Coe * np.exp(-Mas * it),
                               np.arange(itmin, itmax+1),
                               iCorr[i, itmin:itmax+1], sigma=iCorr_err[itmin:itmax+1],
                               p0=np.array((0.1, 0.01)))[0] for i in range(Nconf)])

def fit_effmass_csh (iCorr, itmin, itmax):
    from scipy.optimize import curve_fit
    Nconf     = len(iCorr[:,0])
    iCorr_err = np.std(iCorr, axis=0) * np.sqrt(Nconf-1)
    return np.array([curve_fit(lambda it, Coe, Mas: Coe * np.cosh(Mas * it),
                               np.arange(itmin, itmax+1),
                               iCorr[i, itmin:itmax+1], sigma=iCorr_err[itmin:itmax+1],
                               p0=np.array((0.1, 0.01)))[0] for i in range(Nconf)])

def calc_Zfac (iCorrPS, iCorrSS, imass):
    Nconf = len(iCorrPS[:,0])
    Ntime = len(iCorrPS[0,:])
    ZfacP = np.empty((Nconf, Ntime))
    ZfacS = np.empty((Nconf, Ntime))
    for it in range(Ntime):
        factor = np.sqrt(2.0*imass / np.exp(-imass*it))
        if (np.all(iCorrSS[:,it] > 0)):
            ZfacP[:,it] = iCorrPS[:,it]/np.sqrt(iCorrSS[:,it]) * factor
            ZfacS[:,it] =               np.sqrt(iCorrSS[:,it]) * factor
        else:
            ZfacP[:,it] = np.zeros(Nconf)
            ZfacS[:,it] = np.zeros(Nconf)
    return (ZfacP, ZfacS)

def fit_Zfac_const (iZfac, itmin, itmax):
    Nconf     = len(iZfac[:,0])
    iZfac_err = np.std(iZfac, axis=0) * np.sqrt(Nconf-1)
    return np.array([np.sum(iZfac[i, itmin:itmax+1]/iZfac_err[itmin:itmax+1]) /
                     np.sum(                    1.0/iZfac_err[itmin:itmax+1])
                     for i in range(Nconf)])

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc < 2):
        exit("usage: python %s [ifile(s)1] [ifile(s)2] ..." % basename(argv[0]))
    
    if (main(argv[1:]) != 0):
        exit("ERROR EXIT.")
