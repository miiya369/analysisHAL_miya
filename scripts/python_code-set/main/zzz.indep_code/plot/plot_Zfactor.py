#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc, patches

xrange = [0, 14]
#xrange =  None

yrange = [[-1, 4], [-10e-0, 10e-0], [-0.00015, 0.00005], [0.001, 0.006]]
#yrange = [None,None,None,None]

BinSize     = 1
read_column = 1

tmin = 9
tmax = 11

#lat_spacing = 0.1215
lat_spacing = None

fit_type = "exp"
#fit_type = "cosh"

label_location = 'lower left'

### =========================== Main =========================== ###

def main(ibase, had_name, ofname):
    tmp1 = main_read(ibase, had_name, read_column, BinSize)
    if (tmp1 == -1):
        return -1
    
    tmp2 = main_fit (tmp1[0], tmp1[1], lat_spacing, fit_type, tmin, tmax)
    if (tmp2 == -1):
        return -1
    
    main_calcZ(tmp1[0], tmp1[1], tmp2[2], tmp2[3], tmin, tmax)
    main_plot (tmp1[0], tmp1[1], tmp2[0], tmp2[1], tmp2[2], tmp2[3], tmin, tmax, 
               xrange, yrange, label_location, ofname)
    
    return 0

def main_read(ibase, had_name, rColm, Bsize):
    Corr_PS = read_correlators("%s/correlator.PS.dir/*/%s*" % (ibase, had_name))
    Corr_SS = read_correlators("%s/correlator.SS.dir/*/%s*" % (ibase, had_name))
    Nfile = len(Corr_PS[:,0,0])
    Ntime = len(Corr_PS[0,:,0])
    Ncolm = len(Corr_PS[0,0,:])
    print("# Read correlators: Nfile=%d, Ntime=%d, Ncolm=%d" % (Nfile, Ntime, Ncolm))
    
    if   (Bsize == 0):
        Nbin       = Nfile; 
        Corr_PS_jk = Corr_PS[:,:,rColm]
        Corr_SS_jk = Corr_SS[:,:,rColm]
    elif (Nfile % Bsize != 0 or Bsize == Nfile):
        print("\nERROR: Unexpected binsize.\n"); return -1
    else:
        Nbin       = Nfile // Bsize
        Corr_PS_jk = make_jk_samples(Corr_PS[:,:,rColm], Nbin)
        if (Corr_PS_jk is None):
            return -1
        Corr_SS_jk = make_jk_samples(Corr_SS[:,:,rColm], Nbin)
        if (Corr_SS_jk is None):
            return -1
        print("# Make jackknife samples: Nbin=%d" % Nbin)
    
    return (Corr_PS_jk, Corr_SS_jk)

def main_fit(Corr_PS_jk, Corr_SS_jk, lat_a, ftype, t_min, t_max):
    if (lat_a is None):
        factor_MeV = 1.0
        unit_str   = "[Lattice Unit]"
    else:
        factor_MeV = 197.327053 / lat_a
        unit_str   = "[MeV]"
    
    if   (ftype == "exp"):
        Effmass_PS = calc_effmass_exp(Corr_PS_jk) * factor_MeV
        Effmass_SS = calc_effmass_exp(Corr_SS_jk) * factor_MeV
        fitmass    =  fit_effmass_exp(Corr_PS_jk, t_min, t_max)
    elif (ftype == "cosh"):
        Effmass_PS = calc_effmass_csh(Corr_PS_jk) * factor_MeV
        Effmass_SS = calc_effmass_csh(Corr_SS_jk) * factor_MeV
        fitmass    =  fit_effmass_csh(Corr_PS_jk, t_min, t_max)
    else:
        print("\nERROR: Unknown fit type, %s.\n" % ftype); return -1
    
    #print_data_err(Effmass)
    
    return (Effmass_PS, Effmass_SS, fitmass, (factor_MeV, unit_str))

def main_calcZ(Corr_PS_jk, Corr_SS_jk, fitmass, factor_MeV, t_min, t_max):
    Nbin        =     len(fitmass[:,1])
    fitmass_ave = np.mean(fitmass[:,1])
    fitmass_err = np.std (fitmass[:,1] * factor_MeV[0]) * np.sqrt(Nbin-1)
    
    ZfacP, ZfacS = calc_Zfac(Corr_PS_jk, Corr_SS_jk, fitmass_ave)
    ZfacP       *= factor_MeV[0]
    ZfacS       *= factor_MeV[0]
    fitZfacP     = fit_Zfac_const(ZfacP, t_min, t_max)
    fitZfacS     = fit_Zfac_const(ZfacS, t_min, t_max)
    fitZfacP_ave = np.mean(fitZfacP)
    fitZfacP_err = np.std (fitZfacP) * np.sqrt(Nbin-1)
    fitZfacS_ave = np.mean(fitZfacS)
    fitZfacS_err = np.std (fitZfacS) * np.sqrt(Nbin-1)
    
    fitmass_ave *= factor_MeV[0]
    print("# Effective mass   = %15.16f +/- %15.16f %s" % ( fitmass_ave,  fitmass_err, factor_MeV[1]))
    print("# Z-factor (Point) = %15.16f +/- %15.16f %s" % (fitZfacP_ave, fitZfacP_err, factor_MeV[1]))
    print("# Z-factor (Smear) = %15.16f +/- %15.16f %s" % (fitZfacS_ave, fitZfacS_err, factor_MeV[1]))

def main_plot(Corr_PS_jk, Corr_SS_jk, Effmass_PS, Effmass_SS, fitmass, factor_MeV, t_min, t_max, 
              x_range, y_range, label_loc, ofname):
    rc('text' , usetex=True)
    rc('font' ,**{'family'    : 'Times New Roman',
                  'weight'    : 'bold',
                  'size'      : 12})
    rc('xtick',**{'labelsize' : 14})
    rc('ytick',**{'labelsize' : 14})
    rc('axes' ,**{'labelsize' : 14})
    rc('axes' ,**{'linewidth' : 2})
    
    FigSizeBase = [12, 9]
    FigSize     = [[0.08, 0.58, 0.4, 0.4], [0.58, 0.58, 0.4, 0.4],
                   [0.08, 0.08, 0.4, 0.4], [0.58, 0.08, 0.4, 0.4]]
    ylabel      = ('Effective mass %s','Correlator %s',
                   'Effective Z-factor (Point) %s', 'Effective Z-factor (Smear) %s')
    xlabel      = r'time slice $t-t_0$'
    
    plt.ion(); fig = plt.figure(figsize=FigSizeBase)
    fig.patch.set_facecolor('white')
    axe = [None for i in range(4)]
    for i in range(4):
        axe[i] = fig.add_axes(FigSize[i])
        axe[i].patch.set_facecolor('white')
        axe[i].set_xlabel(xlabel   )
        axe[i].xaxis.set_label_coords( 0.5 , -0.12)
        axe[i].set_ylabel(ylabel[i] % factor_MeV[1])
        axe[i].yaxis.set_label_coords(-0.12,  0.5 )
    
    Ntime       =     len(Corr_PS_jk[0,:])
    Nbin        =     len(fitmass[:,1])
    fitmass_ave = np.mean(fitmass[:,1])
    fitmass_err = np.std (fitmass[:,1] * factor_MeV[0]) * np.sqrt(Nbin-1)

    ZfacP, ZfacS = calc_Zfac(Corr_PS_jk, Corr_SS_jk, fitmass_ave)
    ZfacP       *= factor_MeV[0]
    ZfacS       *= factor_MeV[0]
    fitZfacP     = fit_Zfac_const(ZfacP, t_min, t_max)
    fitZfacS     = fit_Zfac_const(ZfacS, t_min, t_max)
    fitZfacP_ave = np.mean(fitZfacP)
    fitZfacP_err = np.std (fitZfacP) * np.sqrt(Nbin-1)
    fitZfacS_ave = np.mean(fitZfacS)
    fitZfacS_err = np.std (fitZfacS) * np.sqrt(Nbin-1)

    fitmass_ave *= factor_MeV[0]
    
    axe[0].errorbar(np.arange(len(Effmass_PS[0,:]))-0.05, np.mean(Effmass_PS, axis=0), 
                    yerr=(np.std(Effmass_PS, axis=0) * np.sqrt(Nbin-1)), fmt='none', #fmt=None,
                    marker=None, label='Effective mass (PS)', lw=3, ecolor='red', capthick=3, capsize=3)
    axe[0].errorbar(np.arange(len(Effmass_SS[0,:]))+0.05, np.mean(Effmass_SS, axis=0), 
                    yerr=(np.std(Effmass_SS, axis=0) * np.sqrt(Nbin-1)), fmt='none', #fmt=None,
                    marker=None, label='Effective mass (SS)', lw=3, ecolor='blue', capthick=3, capsize=3)
    axe[0].add_patch(patches.Rectangle((t_min, fitmass_ave-fitmass_err), t_max-t_min, 2*fitmass_err, 
                                       color='green', alpha=0.5))
    axe[0].plot((t_min, t_max), (fitmass_ave, fitmass_ave), lw=3, color='green', 
                label='Fit value (t = %d-%d): %6f $\pm$ %6f %s' %
                (t_min, t_max, fitmass_ave, fitmass_err, factor_MeV[1]))
    
    axe[1].errorbar(np.arange(len(Corr_PS_jk[0,:])), np.mean(Corr_PS_jk * factor_MeV[0], axis=0), 
                    yerr=(np.std(Corr_PS_jk * factor_MeV[0], axis=0) * np.sqrt(Nbin-1)), fmt='o', #fmt=None,
                    marker=None, label='Correlator data (PS)', lw=3, color='red', ecolor='red', capthick=3, capsize=3)
    axe[1].errorbar(np.arange(len(Corr_SS_jk[0,:])), np.mean(Corr_SS_jk * factor_MeV[0], axis=0), 
                    yerr=(np.std(Corr_SS_jk * factor_MeV[0], axis=0) * np.sqrt(Nbin-1)), fmt='o', #fmt=None,
                    marker=None, label='Correlator data (SS)', lw=3, color='blue', ecolor='blue', capthick=3, capsize=3)
    t_plot = np.arange(0, Ntime, 0.1)
    axe[1].plot(t_plot, 
                corr_func(t_plot / factor_MeV[0], fitZfacP_ave, fitZfacS_ave, fitmass_ave), 
                lw=3, color='red', label='Fit value (PS): t = %d-%d' % (t_min, t_max))
    axe[1].plot(t_plot, 
                corr_func(t_plot / factor_MeV[0], fitZfacS_ave, fitZfacS_ave, fitmass_ave), 
                lw=3, color='blue', label='Fit value (SS): t = %d-%d' % (t_min, t_max))
    axe[1].set_yscale('log')

    axe[2].errorbar(np.arange(len(ZfacP[0,:])), np.mean(ZfacP, axis=0),
                    yerr=(np.std(ZfacP, axis=0) * np.sqrt(Nbin-1)), fmt='none', #fmt=None,
                    marker=None, label=r'$\sqrt{Z_B^{(\mathrm{Point})}(t)}$', lw=3, ecolor='red', capthick=3, capsize=3)
    axe[2].add_patch(patches.Rectangle((t_min, fitZfacP_ave-fitZfacP_err), t_max-t_min, 2*fitZfacP_err,
                                       color='green', alpha=0.5))
    axe[2].plot((t_min, t_max), (fitZfacP_ave, fitZfacP_ave), lw=3, color='green', 
                label='Fit value (t = %d-%d): %6f $\pm$ %6f %s' %
                (t_min, t_max, fitZfacP_ave, fitZfacP_err, factor_MeV[1]))
    
    axe[3].errorbar(np.arange(len(ZfacS[0,:])), np.mean(ZfacS, axis=0),
                    yerr=(np.std(ZfacS, axis=0) * np.sqrt(Nbin-1)), fmt='none', #fmt=None,
                    marker=None, label=r'$\sqrt{Z_B^{\mathrm{(Smear)}}(t)}$', lw=3, ecolor='blue', capthick=3, capsize=3)
    axe[3].add_patch(patches.Rectangle((t_min, fitZfacS_ave-fitZfacS_err), t_max-t_min, 2*fitZfacS_err,
                                       color='green', alpha=0.5))
    axe[3].plot((t_min, t_max), (fitZfacS_ave, fitZfacS_ave), lw=3, color='green', 
                label='Fit value (t = %d-%d): %6f $\pm$ %6f %s' %
                (t_min, t_max, fitZfacS_ave, fitZfacS_err, factor_MeV[1]))
    
    for i in range(4):
        if (x_range    is not None):
            axe[i].set_xlim(x_range   )
        if (y_range[i] is not None):
            axe[i].set_ylim(y_range[i])
        
        axe[i].grid(which='major',color='gray',linestyle='--')
        #axe[i].axhline(0.0, lw=2, color='black')
        axe[i].legend(numpoints=1, loc=label_loc) #, bbox_to_anchor=(0.0, 0.0))
    
    if (ofname is not None):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    else:
        plt.show(); print("# push Enter to end"); input('# ') #raw_input('# ') # For v.2

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

def corr_func (it, iZ1, iZ2, imass):
    return iZ1 * iZ2 * np.exp(-imass*it) / (2.0 * imass)

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 3 and argc != 4):
        exit("usage: python %s [Path to results] [hadron name] {[ofname]}" % basename(argv[0]))

    ibase    = argv[1].strip()
    had_name = argv[2].strip()
    if (argc == 4):
        ofname = argv[3].strip()
    else:
        ofname = None
        
    if (main(ibase, had_name, ofname) != 0):
        exit("ERROR EXIT.")
