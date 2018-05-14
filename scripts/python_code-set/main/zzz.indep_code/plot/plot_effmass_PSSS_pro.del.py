#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc, patches
rc('text' , usetex=True)
rc('font' ,**{'family'    : 'Times New Roman',
              'weight'    : 'bold',
              'size'      : 14})
rc('xtick',**{'labelsize' : 18})
rc('ytick',**{'labelsize' : 18})
rc('axes' ,**{'labelsize' : 18})
rc('axes' ,**{'linewidth' : 2})
xlim_all = (2, 14)
ylim_all = (0.001, 0.01)

BinSize     = 1
read_column = 1

tmin = 9
tmax = 11

lat_space = 0.1215

### =========================== Main =========================== ###

def main (ibase, ofname):
    Corr_PS_p = read_correlators("%s/correlator.PS.dir/*/proton_CG05_CG05*" % ibase)
    Corr_SS_p = read_correlators("%s/correlator.SS.dir/*/proton_CG05_CG05*" % ibase)
    #Corr_PS_p = read_correlators("%s/correlator.PS.dir/*/proton_NR_NR*" % ibase)
    #Corr_SS_p = read_correlators("%s/correlator.SS.dir/*/proton_NR_NR*" % ibase)
    Corr_PS_d = read_correlators("%s/correlator.PS.dir/*/Delta*" % ibase)
    Corr_SS_d = read_correlators("%s/correlator.SS.dir/*/Delta*" % ibase)
    Nfile = len(Corr_PS_p[:,0,0])
    Ntime = len(Corr_PS_p[0,:,0])
    Ncolm = len(Corr_PS_p[0,0,:])
    print("# Read correlators: Nfile=%d, Ntime=%d, Ncolm=%d" % (Nfile, Ntime, Ncolm))
    
    if   (BinSize == 0):
        Nbin       = Nfile; 
        Corr_PS_p_jk = Corr_PS_p[:,:,read_column]
        Corr_SS_p_jk = Corr_SS_p[:,:,read_column]
        Corr_PS_d_jk = Corr_PS_d[:,:,read_column]
        Corr_SS_d_jk = Corr_SS_d[:,:,read_column]
    elif (Nfile % BinSize != 0 or BinSize == Nfile):
        print("\nERROR: Unexpected binsize.\n"); return -1
    else:
        Nbin       = Nfile // BinSize
        Corr_PS_p_jk = make_jk_samples(Corr_PS_p[:,:,read_column], Nbin)
        if (Corr_PS_p_jk is None):
            return -1
        Corr_PS_d_jk = make_jk_samples(Corr_PS_d[:,:,read_column], Nbin)
        if (Corr_PS_d_jk is None):
            return -1
        Corr_SS_p_jk = make_jk_samples(Corr_SS_p[:,:,read_column], Nbin)
        if (Corr_SS_p_jk is None):
            return -1
        Corr_SS_d_jk = make_jk_samples(Corr_SS_d[:,:,read_column], Nbin)
        if (Corr_SS_d_jk is None):
            return -1
        print("# Make jackknife samples: Nbin=%d" % Nbin)
    
    #Effmass_PS = calc_effmass_exp(Corr_PS_jk)
    #Effmass_SS = calc_effmass_exp(Corr_SS_jk)
    #Effmass_PS = calc_effmass_csh(Corr_PS_jk)
    #Effmass_SS = calc_effmass_csh(Corr_SS_jk)
    
    fitmass_p     = fit_effmass_exp(Corr_PS_p_jk, tmin, tmax)
    fitmass_d     = fit_effmass_exp(Corr_PS_d_jk, tmin, tmax)
    #fitmass     = fit_effmass_csh(Corr_PS_jk, tmin, tmax)
    fitmass_p_ave = np.mean(fitmass_p[:,1])
    fitmass_d_ave = np.mean(fitmass_d[:,1])
    #fitmass_err = np.std (fitmass[:,1]) * np.sqrt(Nbin-1)
    #print("# Effective mass   = %15.16f +/- %15.16f" % (fitmass_ave, fitmass_err))
    
    ZfacP_p, ZfacS_p = calc_Zfac(Corr_PS_p_jk, Corr_SS_p_jk, fitmass_p_ave)
    fitZfacP_p     = fit_Zfac_const(ZfacP_p, tmin, tmax)
    fitZfacS_p     = fit_Zfac_const(ZfacS_p, tmin, tmax)
    fitZfacP_p_ave = np.mean(fitZfacP_p)
    fitZfacP_p_err = np.std (fitZfacP_p) * np.sqrt(Nbin-1)
    fitZfacS_p_ave = np.mean(fitZfacS_p)
    fitZfacS_p_err = np.std (fitZfacS_p) * np.sqrt(Nbin-1)
    
    ZfacP_d, ZfacS_d = calc_Zfac(Corr_PS_d_jk, Corr_SS_d_jk, fitmass_d_ave)
    fitZfacP_d     = fit_Zfac_const(ZfacP_d, tmin, tmax)
    fitZfacS_d     = fit_Zfac_const(ZfacS_d, tmin, tmax)
    fitZfacP_d_ave = np.mean(fitZfacP_d)
    fitZfacP_d_err = np.std (fitZfacP_d) * np.sqrt(Nbin-1)
    fitZfacS_d_ave = np.mean(fitZfacS_d)
    fitZfacS_d_err = np.std (fitZfacS_d) * np.sqrt(Nbin-1)
    
    print("# Proton Z-factor (Point) = %15.16f +/- %15.16f" % (fitZfacP_p_ave, fitZfacP_p_err))
    print("# Proton Z-factor (Smear) = %15.16f +/- %15.16f" % (fitZfacS_p_ave, fitZfacS_p_err))
    print("# Delta  Z-factor (Point) = %15.16f +/- %15.16f" % (fitZfacP_d_ave, fitZfacP_d_err))
    print("# Delta  Z-factor (Smear) = %15.16f +/- %15.16f" % (fitZfacS_d_ave, fitZfacS_d_err))
    
    # --- plot --- #
    plt.ion(); fig = plt.figure()
    fig.patch.set_facecolor('white')
    axe = fig.add_axes((0.15, 0.15, 0.8, 0.75))
    axe.patch.set_facecolor('white')
    
    axe.errorbar(np.arange(len(ZfacP_p[0,:])), np.mean(ZfacP_p, axis=0), 
                 yerr=(np.std(ZfacP_p, axis=0) * np.sqrt(Nbin-1)), fmt=None,
                 marker=None, label=r'$\sqrt{Z_p^{(\mathrm{Point})}}$', lw=3, ecolor='red', capthick=3, capsize=5)
    axe.errorbar(np.arange(len(ZfacP_d[0,:])), np.mean(ZfacP_d, axis=0), 
                 yerr=(np.std(ZfacP_d, axis=0) * np.sqrt(Nbin-1)), fmt=None,
                 marker=None, label=r'$\sqrt{Z_\Delta^{(\mathrm{Point})}}$', lw=3, ecolor='blue', capthick=3, capsize=5)
    
    #axe.errorbar(np.arange(len(ZfacS_p[0,:])), np.mean(ZfacS_p, axis=0), 
    #             yerr=(np.std(ZfacS_p, axis=0) * np.sqrt(Nbin-1)), fmt=None,
    #             marker=None, label=r'$\sqrt{Z_p^{(\mathrm{Smear})}}$', lw=3, ecolor='red', capthick=3, capsize=5)
    #axe.errorbar(np.arange(len(ZfacS_d[0,:])), np.mean(ZfacS_d, axis=0), 
    #             yerr=(np.std(ZfacS_d, axis=0) * np.sqrt(Nbin-1)), fmt=None,
    #             marker=None, label=r'$\sqrt{Z_\Delta^{(\mathrm{Smear})}}$', lw=3, ecolor='blue', capthick=3, capsize=5)
    
    axe.add_patch(patches.Rectangle((tmin, fitZfacP_p_ave-fitZfacP_p_err), tmax-tmin, 2*fitZfacP_p_err, 
                                     color='red', alpha=0.5))
    axe.plot((tmin, tmax), (fitZfacP_p_ave, fitZfacP_p_ave), lw=3, color='red', label='%1.5f'%fitZfacP_p_ave)
    axe.add_patch(patches.Rectangle((tmin, fitZfacP_d_ave-fitZfacP_d_err), tmax-tmin, 2*fitZfacP_d_err, 
                                     color='blue', alpha=0.5))
    axe.plot((tmin, tmax), (fitZfacP_d_ave, fitZfacP_d_ave), lw=3, color='blue', label='%1.5f'%fitZfacP_d_ave)
    
    #axe.add_patch(patches.Rectangle((tmin, fitZfacS_p_ave-fitZfacS_p_err), tmax-tmin, 2*fitZfacS_p_err, 
    #                                 color='red', alpha=0.5))
    #axe.plot((tmin, tmax), (fitZfacS_p_ave, fitZfacS_p_ave), lw=3, color='red', label='%1.5f'%fitZfacS_p_ave)
    #axe.add_patch(patches.Rectangle((tmin, fitZfacS_d_ave-fitZfacS_d_err), tmax-tmin, 2*fitZfacS_d_err, 
    #                                 color='blue', alpha=0.5))
    #axe.plot((tmin, tmax), (fitZfacS_d_ave, fitZfacS_d_ave), lw=3, color='blue', label='%1.5f'%fitZfacS_d_ave)
    
    if (ylim_all is not None):
        axe.set_ylim(ylim_all)
    if (xlim_all is not None):
        axe.set_xlim(xlim_all)
    #axe.set_title(r'Fit range = [%d,%d], Mass = %1.5f * (197.327/%1.4f) = %4.2f MeV' % 
    #              (tmin, tmax, fitmass_ave, lat_space, (fitmass_ave*197.327/lat_space)), x=0.45, y=1.05)
    axe.legend(numpoints=1, loc='upper right')
    axe.set_ylabel(r'$\sqrt{Z_B^{(\mathrm{Smear})}}$')
    axe.set_xlabel(r'time slice ($t-t_0$)')
    axe.yaxis.set_label_coords(-0.1,0.5)
    axe.xaxis.set_label_coords(0.5,-0.1)
    #axe.axhline(0.0, lw=2, color='black')
    axe.grid(which='major',color='gray',linestyle='--')
    plt.show(); print("# push Enter to end"); raw_input('# ') # For v.2
    if (ofname is not None):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    
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
    return np.array([curve_fit(lambda it, Coe, Mas: Coe * np.cosh(-Mas * it),
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
                     np.sum(                    1.0/iZfac_err[itmin:itmax+1]) for i in range(Nconf)])

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 2 and argc != 3):
        exit("usage: python %s [Path to results] {[ofname]}" % basename(argv[0]))

    ibase = argv[1].strip()
    if (argc == 3):
        ofname = argv[2].strip()
    else:
        ofname = None
    
    if (main(ibase, ofname) != 0):
        exit("ERROR EXIT.")
