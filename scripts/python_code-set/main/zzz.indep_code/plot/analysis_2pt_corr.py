#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Takaya Miyamoto (YITP)
# Brief : analysis 2pt correlators
# Date  : Thu May 17 16:19:40 JST 2018

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc, patches

### ================== Global Parameters Init. ================= ###
rc('text' , usetex=True)
rc('font' ,**{'family'    : 'Times New Roman',
              'weight'    : 'bold',
              'size'      : 11})
rc('xtick',**{'labelsize' : 12})
rc('ytick',**{'labelsize' : 12})
rc('axes' ,**{'labelsize' : 12})
rc('axes' ,**{'linewidth' : 2})

read_column = 1
lattice_spacing = None

xlim_common  = (5, 25)
ylim_effmass = None
ylim_corr    = None
ylim_ZfacP   = None
ylim_ZfacS   = None
#ylim_effmass = (0.5, 1.0)
#ylim_corr    = (1e-16, 1e-2)
#ylim_ZfacP   = (0.02, 0.1)
#ylim_ZfacS   = (0.00002, 0.0001)

label_location = 'upper left'
### ============================================================ ###



###==============================================###
###=== Choose main part which you want to use ===###
def main(argc, argv): 
    #return main_print_correlator(argc, argv)
    #return main_print_effmass(argc, argv)
    #return main_print_Zfactor(argc, argv)
    #return main_fit_effmass(argc, argv)
    #return main_fit_Zfactor(argc, argv)
    #return main_plot_effmass(argc, argv)
    #return main_plot_effmass_PS_SS(argc, argv)
    #return main_plot_Zfactor(argc, argv)
    print("Warning: Choose main part which you want to use."); return -1
###==============================================###



### ============================================================ ###
### =========================== Main =========================== ###
def main_print_correlator(argc, argv):
    from os.path import basename
    if (argc < 3):
        print("usage: python", basename(argv[0]), 
              "[binsize] [ifile(s)1] [ifile(s)2] ...")
        return 0
    bsize = int(argv[1])
    ipath =     argv[2:]
    print_data_err(read_correlators_jk(ipath, bsize, read_column))
    return 0

def main_print_effmass(argc, argv):
    from os.path import basename
    if (argc < 4):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [ifile(s)1] [ifile(s)2] ...")
        return 0
    ftype =     argv[1]
    bsize = int(argv[2])
    ipath =     argv[3:]
    
    if (lattice_spacing is None):
        factor_MeV = 1.0
    else:
        factor_MeV = 197.327053 / lattice_spacing
    
    print_data_err(calc_effmass(read_correlators_jk(ipath, bsize, read_column), ftype) * factor_MeV)
    return 0

def main_print_Zfactor(argc, argv):
    from os.path import basename
    if (argc != 7):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [tmin_fit] [tmax_fit] [Path to results] [hadron name]")
        return 0
    ftype    =     argv[1]
    bsize    = int(argv[2])
    tmin_fit = int(argv[3])
    tmax_fit = int(argv[4])
    ibase    =     argv[5]
    had_name =     argv[6]
    Zfacs = calc_Zfactor(*read_correlators_PS_SS_jk(ibase, had_name, bsize, read_column), tmin_fit, tmax_fit, ftype)
    print(    "# Zfac_P:"); print_data_err(Zfacs[0,:,:])
    print("\n\n# Zfac_S:"); print_data_err(Zfacs[1,:,:])
    return 0

def main_fit_effmass(argc, argv):
    from os.path import basename
    if (argc < 6):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [tmin_fit] [tmax_fit] [ifile(s)1] [ifile(s)2] ...")
        return 0
    ftype    =     argv[1]
    bsize    = int(argv[2])
    tmin_fit = int(argv[3])
    tmax_fit = int(argv[4])
    ipath    =     argv[5:]
    fit_effmass(read_correlators_jk(ipath, bsize, read_column), tmin_fit, tmax_fit, lattice_spacing, ftype)
    return 0

def main_fit_Zfactor(argc, argv):
    from os.path import basename
    if (argc != 9):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [tmin_fit_mass] [tmax_fit_mass] [tmin_fit_Zfac] [tmax_fit_Zfac] "+
              "[Path to results] [hadron name]")
        return 0
    ftype         =     argv[1]
    bsize         = int(argv[2])
    tmin_fit_mass = int(argv[3])
    tmax_fit_mass = int(argv[4])
    tmin_fit_Zfac = int(argv[5])
    tmax_fit_Zfac = int(argv[6])
    ibase         =     argv[7]
    had_name      =     argv[8]
    fit_Zfactor(*read_correlators_PS_SS_jk(ibase, had_name, bsize, read_column), 
                 tmin_fit_mass, tmax_fit_mass, tmin_fit_Zfac, tmax_fit_Zfac, ftype)
    return 0

def main_plot_effmass(argc, argv):
    from os.path import basename
    if (argc < 7):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [tmin_fit] [tmax_fit] [ofname (or 'none')] [ifile(s)1] [ifile(s)2] ...")
        return 0
    ftype    =     argv[1]
    bsize    = int(argv[2])
    tmin_fit = int(argv[3])
    tmax_fit = int(argv[4])
    ofname   =     argv[5]
    ipath    =     argv[6:]
    
    if (lattice_spacing is None):
        factor_MeV = 1.0
        unit_str   = "Lattice Unit"
    else:
        factor_MeV = 197.327053 / lattice_spacing
        unit_str   = "MeV"
    
    icorr   = read_correlators_jk(ipath, bsize, read_column)
    effmass = calc_effmass(icorr, ftype) * factor_MeV
    fitmass =  fit_effmass(icorr, tmin_fit, tmax_fit, lattice_spacing, ftype)
    
    plt.ion(); fig = plt.figure()
    fig.patch.set_facecolor('white')
    axe = fig.add_axes((0.15, 0.15, 0.8, 0.8))
    plot_data_with_fit(axe, effmass, fitmass, tmin_fit, tmax_fit, 
                       xlim_common, ylim_effmass, r"Effective mass [%s]" % unit_str)
    if (ofname != "none"):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    else:
        plt.show(); print("# push Enter to end"); input('# ') #raw_input('# ') # For v.2
    return 0

def main_plot_effmass_PS_SS(argc, argv):
    from os.path import basename
    if (argc != 8):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [tmin_fit] [tmax_fit] [Path to results] [hadron name] [ofname (or 'none')]")
        return 0
    ftype    =     argv[1]
    bsize    = int(argv[2])
    tmin_fit = int(argv[3])
    tmax_fit = int(argv[4])
    ibase    =     argv[5]
    had_name =     argv[6]
    ofname   =     argv[7]
    
    if (lattice_spacing is None):
        factor_MeV = 1.0
        unit_str   = "Lattice Unit"
    else:
        factor_MeV = 197.327053 / lattice_spacing
        unit_str   = "MeV"
    
    icorr_PS, icorr_SS = read_correlators_PS_SS_jk(ibase, had_name, bsize, read_column)
    effmass_PS = calc_effmass(icorr_PS, ftype) * factor_MeV
    effmass_SS = calc_effmass(icorr_SS, ftype) * factor_MeV
    fitmass_PS =  fit_effmass(icorr_PS, tmin_fit, tmax_fit, lattice_spacing, ftype)
    
    plt.ion(); fig = plt.figure()
    fig.patch.set_facecolor('white')
    axe = fig.add_axes((0.15, 0.15, 0.8, 0.8))
    plot_effmass_PS_SS(axe, effmass_PS, effmass_SS, fitmass_PS, tmin_fit, tmax_fit,
                       xlim_common, ylim_effmass, unit_str)
    if (ofname != "none"):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    else:
        plt.show(); print("# push Enter to end"); input('# ') #raw_input('# ') # For v.2
    return 0

def main_plot_Zfactor(argc, argv):
    from os.path import basename
    if (argc != 10):
        print("usage: python", basename(argv[0]), 
              "['exp' or 'cosh'] [binsize] [tmin_fit_mass] [tmax_fit_mass] [tmin_fit_Zfac] [tmax_fit_Zfac] "+
              "[Path to results] [hadron name] [ofname (or 'none')]")
        return 0
    ftype         =     argv[1]
    bsize         = int(argv[2])
    tmin_fit_mass = int(argv[3])
    tmax_fit_mass = int(argv[4])
    tmin_fit_Zfac = int(argv[5])
    tmax_fit_Zfac = int(argv[6])
    ibase         =     argv[7]
    had_name      =     argv[8]
    ofname        =     argv[9]
    
    if (lattice_spacing is None):
        factor_MeV = 1.0
        unit_str   = "Lattice Unit"
    else:
        factor_MeV = 197.327053 / lattice_spacing
        unit_str   = "MeV"
    
    icorr_PS, icorr_SS = read_correlators_PS_SS_jk(ibase, had_name, bsize, read_column)
    effmass_PS = calc_effmass(icorr_PS, ftype) * factor_MeV
    effmass_SS = calc_effmass(icorr_SS, ftype) * factor_MeV
    fitmass_PS =  fit_effmass(icorr_PS, tmin_fit_mass, tmax_fit_mass, lattice_spacing, ftype)
    Zfacs      = calc_Zfactor(icorr_PS, icorr_SS, tmin_fit_mass, tmax_fit_mass, ftype)
    fitZfacs   =  fit_Zfactor(icorr_PS, icorr_SS, tmin_fit_mass, tmax_fit_mass, 
                              tmin_fit_Zfac, tmax_fit_Zfac, ftype)
    
    plt.ion(); fig = plt.figure(figsize=(12, 9))
    fig.patch.set_facecolor('white')
    axe = [fig.add_axes((0.08, 0.58, 0.4, 0.4)),
           fig.add_axes((0.58, 0.58, 0.4, 0.4)),
           fig.add_axes((0.08, 0.08, 0.4, 0.4)),
           fig.add_axes((0.58, 0.08, 0.4, 0.4))]
    plot_effmass_PS_SS(axe[0], effmass_PS, effmass_SS, fitmass_PS, tmin_fit_mass, tmax_fit_mass,
                       xlim_common, ylim_effmass, unit_str)
    plot_2ptcorr_with_fit(axe[1], icorr_PS, icorr_SS, np.mean(fitmass_PS/factor_MeV),
                          np.mean(fitZfacs[0]), np.mean(fitZfacs[1]), xlim_common, ylim_corr, ftype)
    axe[1].set_yscale('log')
    plot_data_with_fit(axe[2], Zfacs[0], fitZfacs[0], tmin_fit_Zfac, tmax_fit_Zfac, 
                       xlim_common, ylim_ZfacP, r"Z-factor (Point)")
    plot_data_with_fit(axe[3], Zfacs[1], fitZfacs[1], tmin_fit_Zfac, tmax_fit_Zfac, 
                       xlim_common, ylim_ZfacS, r"Z-factor (Smear)")
    if (ofname != "none"):
        print("# Output to '%s'" % ofname); plt.savefig(ofname)
    else:
        plt.show(); print("# push Enter to end"); input('# ') #raw_input('# ') # For v.2
    return 0
### ======================== Functions ========================= ###
def read_correlator_from_path(ifpath, auto_fb = True):
    from glob    import glob, iglob
    from os.path import dirname, basename
    icorr       = np.array([np.loadtxt(ifname) for ifname in iglob(ifpath)])
    ifpath_anti = dirname(ifpath) + "/anti" + basename(ifpath)
    if (len(glob(ifpath_anti)) == 0 or not auto_fb):
        return icorr
    else:
        print("# Read anti-hadron...")
        icorr_anti = np.array([np.loadtxt(ifname) for ifname in iglob(ifpath_anti)])
        return (icorr + np.roll(icorr_anti[:,::-1,:], 1, axis=1)) / 2.0

def read_correlators(ifpaths, auto_fb = True, verbose_flg = True):
    icorrs = read_correlator_from_path(ifpaths[0], auto_fb = auto_fb)
    for i in range(1, len(ifpaths)):
        icorrs = np.append(icorrs, read_correlator_from_path(ifpaths[i], auto_fb=auto_fb), axis=0)
    if (verbose_flg):
        Nfile = len(icorrs[:,0,0])
        Ntime = len(icorrs[0,:,0])
        Ncolm = len(icorrs[0,0,:])
        print("# Read correlators: Nfile=%d, Ntime=%d, Ncolm=%d" % (Nfile, Ntime, Ncolm))
    return icorrs

def make_jk_samples (idata, binsize, verbose_flg = True):
    Ndata = len(idata[:,0])
    if   (binsize == 0):
        return idata
    elif (Ndata % binsize != 0 or binsize == Ndata):
        print("\nERROR: Unexpected binsize, #.data(=%d) with binsize(=%d).\n" % (Ndata, binsize))
        return None
    else:
        Nbin = Ndata // binsize
        if (verbose_flg):
            print("# Make jackknife samples: Nbin=%d" % Nbin)
        return np.array([(np.sum(idata, axis=0)-np.sum(idata[i*binsize:(i+1)*binsize, :], axis=0))
                         for i in range(Nbin)]) / float(Ndata-binsize)

def print_data_err (idata, is_jk_data = True):
    Ndata = len(idata[:,0])
    if (is_jk_data):
        factor_jk = np.sqrt(Ndata-1)
    else:
        factor_jk = 1.0 / np.sqrt(Ndata-1)
    idata_ave = np.mean(idata, axis=0)
    idata_err = np.std (idata, axis=0) * factor_jk
    for i in range(len(idata_ave)):
        print("%d %15.16e %15.16e" % (i, idata_ave[i], idata_err[i]))

def calc_effmass(icorr, effmass_type = "exp"):
    Nconf   = len(icorr[:,0])
    Ntime   = len(icorr[0,:])
    effmass = np.zeros((Nconf, Ntime-1))
    for it in range(Ntime-1):
        if (np.all(icorr[:,it] != 0)):
            if   (effmass_type == "exp"):
                tmp_d = icorr[:,it+1] / icorr[:,it]
                if (np.all(tmp_d > 0)):
                    effmass[:,it] = -np.log(tmp_d)
            elif (effmass_type == "cosh"):
                tmp_d = (icorr[:,it+1] + icorr[:,it-1]) / (2.0*icorr[:,it])
                if (np.all(tmp_d > 1)):
                    effmass[:,it] = np.arccosh(tmp_d)
            else:
                print("\nERROR: Unknown type, %s.\n" % effmass_type); return None
    return effmass

def read_correlators_jk(ifpaths, binsize, read_colm = 1):
    return make_jk_samples(read_correlators(ifpaths)[:,:,read_column], binsize)

def read_correlators_PS_SS_jk(ibase, had_name, binsize, read_colm = 1):
    return (read_correlators_jk(("%s/correlator.PS.dir/*/%s*" % (ibase, had_name),),
                                binsize, read_colm),
            read_correlators_jk(("%s/correlator.SS.dir/*/%s*" % (ibase, had_name),),
                                binsize, read_colm))

def fit_effmass(icorr, tmin_fit, tmax_fit, lat_spacing = None, effmass_type = "exp",
                is_jk_data = True, verbose_flg = True):
    
    from scipy.optimize import curve_fit
    Nconf = len(icorr[:,0])
    harfT = len(icorr[0,:]) // 2
    if (is_jk_data):
        factor_jk = np.sqrt(Nconf-1)
    else:
        factor_jk = 1.0 / np.sqrt(Nconf-1)
    icorr_err = np.std(icorr, axis=0) * factor_jk
    
    if   (effmass_type == "exp"):
        fit_func = lambda it, coe, mas: coe * np.exp(-mas *  it)
    elif (effmass_type == "cosh"):
        fit_func = lambda it, coe, mas: coe * np.cosh(mas * (it-harfT))
    else:
        print("\nERROR: Unknown type, %s.\n" % effmass_type); return None
    
    if (lat_spacing is None):
        factor_MeV = 1.0
        unit_str   = "Lattice Unit"
    else:
        factor_MeV = 197.327053 / lat_spacing
        unit_str   = "MeV"
    
    fitmass = np.array([curve_fit(fit_func, np.arange(tmin_fit, tmax_fit+1),
                                  icorr[i, tmin_fit:tmax_fit+1], sigma=icorr_err[tmin_fit:tmax_fit+1],
                                  p0=np.array((0.1, 0.01)))[0][1] for i in range(Nconf)]) * factor_MeV
    fitmass = np.fabs(fitmass)
    if (verbose_flg):
        fitmass_ave = np.mean(fitmass)
        fitmass_err = np.std (fitmass) * factor_jk
        print("# Effective mass   = %15.16f +/- %15.16f [%s]" % (fitmass_ave, fitmass_err, unit_str))
    return fitmass

def calc_Zfactor(icorr_PS, icorr_SS, tmin_fit, tmax_fit, effmass_type = "exp"):
    Nconf   = len(icorr_PS[:,0])
    Ntime   = len(icorr_PS[0,:])
    harfT   = Ntime // 2
    fitmass = fit_effmass(icorr_PS, tmin_fit, tmax_fit, None, effmass_type, verbose_flg=False)
    
    if   (effmass_type == "exp"):
        factor = np.array([[2.0*fitmass[ic] /  np.exp(-fitmass[ic]*it) 
                            for it in range(Ntime)] for ic in range(Nconf)])
    elif (effmass_type == "cosh"):
        factor = np.array([[fitmass[ic] / (np.exp(-fitmass[ic]*harfT)*np.cosh(fitmass[ic]*(it-harfT)))
                            for it in range(Ntime)] for ic in range(Nconf)])
    else:
        print("\nERROR: Unknown type, %s.\n" % effmass_type); return None
    
    ZpZs = icorr_PS * factor
    ZsZs = icorr_SS * factor
    
    Zs = np.array([[np.sqrt(ZsZs[ic,it])  if (ZsZs[ic,it]  > 0) else 0.0 for it in range(Ntime)] for ic in range(Nconf)])
    Zp = np.array([[ZpZs[ic,it]/Zs[ic,it] if (  Zs[ic,it] != 0) else 0.0 for it in range(Ntime)] for ic in range(Nconf)])
    
    return np.array([Zp, Zs])

def fit_Zfactor(icorr_PS, icorr_SS, tmin_fit_mass, tmax_fit_mass, tmin_fit_Zfac, tmax_fit_Zfac, 
                effmass_type = "exp", is_jk_data = True, verbose_flg = True):
    
    Zfacs = calc_Zfactor(icorr_PS, icorr_SS, tmin_fit_mass, tmax_fit_mass, effmass_type)
    Nconf = len(Zfacs[0,:,0])
    if (is_jk_data):
        factor_jk = np.sqrt(Nconf-1)
    else:
        factor_jk = 1.0 / np.sqrt(Nconf-1)
    
    Zfacs_err = np.array([ np.std(Zfacs[ips,:,:], axis=0) * factor_jk for ips in range(2) ])
    Zfacs_fit = np.array([[np.sum(Zfacs[ips,ic,tmin_fit_Zfac:tmax_fit_Zfac+1]/Zfacs_err[ips,tmin_fit_Zfac:tmax_fit_Zfac+1]) /
                           np.sum(                                        1.0/Zfacs_err[ips,tmin_fit_Zfac:tmax_fit_Zfac+1])
                           for ic in range(Nconf)] for ips in range(2)])
    if (verbose_flg):
        Zfacs_ave = np.mean(Zfacs_fit, axis=1)
        Zfacs_err = np.std (Zfacs_fit, axis=1) * factor_jk
        print("# Z-factor (Point) = %15.16f +/- %15.16f" % (Zfacs_ave[0], Zfacs_err[0]))
        print("# Z-factor (Smear) = %15.16f +/- %15.16f" % (Zfacs_ave[1], Zfacs_err[1]))
    return Zfacs_fit

def plot_data_with_fit(axe, idata, fitdata, tmin_fit, tmax_fit, xlim_in, ylim_in, label_in, is_jk_data = True):
    Nconf = len(idata[:,0])
    Ntime = len(idata[0,:])
    if (is_jk_data):
        factor_jk = np.sqrt(Nconf-1)
    else:
        factor_jk = 1.0 / np.sqrt(Nconf-1)
    fitdata_ave = np.mean(fitdata)
    fitdata_err = np.std (fitdata) * factor_jk
    
    axe.patch.set_facecolor('white')
    axe.errorbar(np.arange(xlim_in[0],xlim_in[1]+1), np.mean(idata[:,xlim_in[0]:xlim_in[1]+1], axis=0),
                 yerr=(np.std(idata[:,xlim_in[0]:xlim_in[1]+1], axis=0) * factor_jk), fmt='none',
                 marker=None, label=r'Original data', lw=3, ecolor='red', capthick=3, capsize=3)
    axe.add_patch(patches.Rectangle((tmin_fit, fitdata_ave-fitdata_err), tmax_fit-tmin_fit, 2*fitdata_err, 
                                    color='green', alpha=0.5))
    axe.plot((tmin_fit, tmax_fit), (fitdata_ave, fitdata_ave), lw=3, color='green',
             label=r'Fit value ($t = %d-%d$): %6f $\pm$ %6f' % (tmin_fit, tmax_fit, fitdata_ave, fitdata_err))
    
    if (ylim_in is not None):
        axe.set_ylim(ylim_in)
    if (xlim_in is not None):
        axe.set_xlim(xlim_in)
    axe.legend(numpoints=1, loc=label_location)
    axe.set_ylabel(label_in)
    axe.set_xlabel(r'time slice ($t-t_0$)')
    axe.yaxis.set_label_coords(-0.12,0.5)
    axe.xaxis.set_label_coords(0.5,-0.1)
    axe.grid(which='major',color='gray',linestyle='--')

def plot_effmass_PS_SS(axe, effmass_PS, effmass_SS, fitmass, tmin_fit, tmax_fit, 
                       xlim_in, ylim_in, unit_str = "Lattice Unit", is_jk_data = True):
    Nconf = len(effmass_PS[:,0])
    Ntime = len(effmass_PS[0,:])
    if (is_jk_data):
        factor_jk = np.sqrt(Nconf-1)
    else:
        factor_jk = 1.0 / np.sqrt(Nconf-1)
    fitmass_ave = np.mean(fitmass)
    fitmass_err = np.std (fitmass) * factor_jk
    
    axe.patch.set_facecolor('white')
    axe.errorbar(np.arange(xlim_in[0],xlim_in[1]+1)-0.1, np.mean(effmass_PS[:,xlim_in[0]:xlim_in[1]+1], axis=0),
                 yerr=(np.std(effmass_PS[:,xlim_in[0]:xlim_in[1]+1], axis=0) * factor_jk), fmt='none',
                 marker=None, label=r'Effective mass data (PS)', lw=3, ecolor='red', capthick=3, capsize=3)
    axe.errorbar(np.arange(xlim_in[0],xlim_in[1]+1)+0.1, np.mean(effmass_SS[:,xlim_in[0]:xlim_in[1]+1], axis=0),
                 yerr=(np.std(effmass_SS[:,xlim_in[0]:xlim_in[1]+1], axis=0) * factor_jk), fmt='none',
                 marker=None, label=r'Effective mass data (SS)', lw=3, ecolor='blue', capthick=3, capsize=3)
    axe.add_patch(patches.Rectangle((tmin_fit, fitmass_ave-fitmass_err), tmax_fit-tmin_fit, 2*fitmass_err, 
                                    color='green', alpha=0.5))
    axe.plot((tmin_fit, tmax_fit), (fitmass_ave, fitmass_ave), lw=3, color='green',
             label=r'Fit value ($t = %d-%d$): %6f $\pm$ %6f [%s]' % (tmin_fit, tmax_fit, fitmass_ave, fitmass_err, unit_str))
    
    if (ylim_in is not None):
        axe.set_ylim(ylim_in)
    if (xlim_in is not None):
        axe.set_xlim(xlim_in)
    axe.legend(numpoints=1, loc=label_location)
    axe.set_ylabel(r'Effective mass [%s]' % unit_str)
    axe.set_xlabel(r'time slice ($t-t_0$)')
    axe.yaxis.set_label_coords(-0.12,0.5)
    axe.xaxis.set_label_coords(0.5,-0.1)
    axe.grid(which='major',color='gray',linestyle='--')

def plot_2ptcorr_with_fit(axe, icorr_PS, icorr_SS, fitmass_ave, fitZfacP_ave, fitZfacS_ave,
                          xlim_in, ylim_in, effmass_type = "exp", is_jk_data = True):
    Nconf = len(icorr_PS[:,0])
    Ntime = len(icorr_PS[0,:])
    harfT = Ntime // 2
    if (is_jk_data):
        factor_jk = np.sqrt(Nconf-1)
    else:
        factor_jk = 1.0 / np.sqrt(Nconf-1)
    
    if   (effmass_type == "exp"):
        corr_func = lambda it, Z1, Z2, m: Z1 * Z2 * np.exp(-m*it) / (2.0*m)
    elif (effmass_type == "cosh"):
        corr_func = lambda it, Z1, Z2, m: Z1 * Z2 * np.exp(-m*harfT) * np.cosh(m*(it-harfT)) / m
    else:
        print("\nERROR: Unknown type, %s.\n" % effmass_type); return None
    
    axe.patch.set_facecolor('white')
    axe.errorbar(np.arange(xlim_in[0],xlim_in[1]+1), np.mean(icorr_PS[:,xlim_in[0]:xlim_in[1]+1], axis=0),
                 yerr=(np.std(icorr_PS[:,xlim_in[0]:xlim_in[1]+1], axis=0) * factor_jk), fmt='ro',
                 marker=None, label=r'Correlator data (PS)', lw=3, ecolor='red', capthick=3, capsize=3)
    axe.errorbar(np.arange(xlim_in[0],xlim_in[1]+1), np.mean(icorr_SS[:,xlim_in[0]:xlim_in[1]+1], axis=0),
                 yerr=(np.std(icorr_SS[:,xlim_in[0]:xlim_in[1]+1], axis=0) * factor_jk), fmt='bo',
                 marker=None, label=r'Correlator data (SS)', lw=3, ecolor='blue', capthick=3, capsize=3)
    tplot = np.arange(xlim_in[0],xlim_in[1]+1, 0.1)
    axe.plot(tplot, corr_func(tplot,fitZfacP_ave,fitZfacS_ave,fitmass_ave), lw=3, color='magenta', label=r'Fit value (PS)')
    axe.plot(tplot, corr_func(tplot,fitZfacS_ave,fitZfacS_ave,fitmass_ave), lw=3, color='cyan', label=r'Fit value (SS)')
    
    if (ylim_in is not None):
        axe.set_ylim(ylim_in)
    if (xlim_in is not None):
        axe.set_xlim(xlim_in)
    axe.legend(numpoints=1, loc=label_location)
    axe.set_ylabel(r'2pt-Correlator')
    axe.set_xlabel(r'time slice ($t-t_0$)')
    axe.yaxis.set_label_coords(-0.12,0.5)
    axe.xaxis.set_label_coords(0.5,-0.1)
    axe.grid(which='major',color='gray',linestyle='--')
### ============================================================ ###
### ============================================================ ###
if (__name__ == "__main__"):
    from sys import exit, argv
    if (main(len(argv), argv) != 0):
        exit("ERROR EXIT.")
