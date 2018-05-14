#!/usr/bin/python

from __future__ import print_function

import numpy as np
from   glob  import glob
import matplotlib.pyplot as plt
from   matplotlib import rc
rc('text' , usetex=True)
rc('font' ,**{'family'    : 'Times New Roman',
              'weight'    : 'bold',
              'size'      : 11})
rc('xtick',**{'labelsize' : 12})
rc('ytick',**{'labelsize' : 12})
rc('axes' ,**{'labelsize' : 12})
rc('axes' ,**{'linewidth' : 2})
xlim_all = None
ylim_all = None
#xlim_all = (0, 14)
#ylim_all = (-6e-24, 1.4e-23)

gLsize  = 16
gTime   = 6
gSpin   = 0
gTshift = "000"
#gTshift = "ALL"
gSprj_flg = True

gNBSch   = ("Prot__Lamb_", "Prot__SigZ_", "Neut__SigP_")
gSrc_op  = "CG05"
gBase    = "/Users/miiya/Dropbox/programs/data/hf_check_old_new/res_new"

ispinU = np.array([
        [[1.0,              0.0,               0.0],
         [0.0, -np.sqrt(1.0/3.0), +np.sqrt(2.0/3.0)],
         [0.0, +np.sqrt(2.0/3.0), +np.sqrt(1.0/3.0)]],
        
        [[1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0],
         [0.0, 0.0, 1.0]]
        ])

flavorU = np.array([
        [[+np.sqrt(9.0/10.0), -np.sqrt(1.0/10.0), 0.0],
         [-np.sqrt(1.0/10.0), -np.sqrt(9.0/10.0), 0.0],
         [              0.0 ,               0.0 , 1.0]],
        
        [[+np.sqrt(1.0/2.0), -np.sqrt(1.0/2.0), 0.0],
         [-np.sqrt(1.0/2.0), -np.sqrt(1.0/2.0), 0.0],
         [             0.0 ,              0.0 , 1.0]]
        ])

### =========================== Main =========================== ###

def main (fconf_lst):
    with open(fconf_lst, 'r') as clst:
        confs = clst.readlines() #; print(confs); return -1
    Nconf = len(confs) ; print("# N.conf = %d" % Nconf)
    Nch   = len(gNBSch)
    idata = np.array([read_NBS_NxN_oct(gBase, gNBSch, gSpin, confs[i].strip(), gTime, gTshift, 
                                       gSrc_op, gSprj_flg, gLsize) for i in range(Nconf)])
    if (Nconf != 1):
        print("# Make Jackknife samples...")
        Ndata = len(idata[0,:,0,0])
        idata = np.reshape(make_jk_samples(np.reshape(idata, (Nconf,Ndata*Nch*Nch)), Nconf), (Nconf,Ndata,Nch,Nch))
    
    idata = np.array([channel_projection(idata[i,:,:,:], ispinU [    0,:,:]) for i in range(Nconf)])
    idata = np.array([channel_projection(idata[i,:,:,:], flavorU[gSpin,:,:]) for i in range(Nconf)])
    
    Rcood_plot, idata_ave_plot, idata_err_plot = set_data_plot(idata, gLsize)
    plt.ion(); fig, axe = plt.subplots(nrows=Nch, ncols=Nch, figsize=(12,8))
    plt.subplots_adjust(wspace=0.2, hspace=0.2); fig.patch.set_facecolor('white')
    plot_data_NxN(axe, Rcood_plot, idata_ave_plot, idata_err_plot, 'NBSwave', 1, '.', 'red')
    #plot_data(axe[0,1], Rcood_plot, idata_ave_plot[:,1,1], idata_err_plot[:,1,1], 'NBSwave', 1, '.', 'blue', xlim_all, ylim_all)
    plt.show(); print("# push Enter to end"); input('# ')
    #plt.savefig('fig.eps')
    
    return 0

### ============================================================ ###

def read_NBS_NxN_oct (ipath_dir, iNBSch, ispin, iconf, itime, itshift, isrc_op, sprj_flg, iLsize):
    Nch = len(iNBSch)
    
    if (sprj_flg):
        ifiles = [[glob("%s/NBS_2Boct.%s_%s.dir/%s/NBSwave.%+04d+%s.000.000.000.%s.*_CG05.*_%s" % 
                        (ipath_dir, iNBSch[i], iNBSch[j], iconf, itime, itshift, iconf, isrc_op)) 
                   for j in range(Nch)] for i in range(Nch)]
    else:
        ifiles = [[glob("%s/Proj.NBS_2Boct.%s_%s.dir/spin%d_ave/%s/NBSwave.%+04d+%s.000.000.000.%s.*_CG05.*_%s" %
                        (ipath_dir, iNBSch[i], iNBSch[j], ispin, iconf, itime, itshift, iconf, isrc_op)) 
                   for j in range(Nch)] for i in range(Nch)]
    
    for i in range(Nch):
        for j in range(Nch):
            if(len(ifiles[i][j]) != 1):
                print("\nERROR: Unexpected files in:"); print(ifiles); return None
    
    idata = np.array([[np.fromfile(ifiles[i][j][0], '>d') for j in range(Nch)] for i in range(Nch)])
    
    if (sprj_flg):
        print("# Spin projection...")
        idata = np.array([[spin_proj_oct(idata[i,j,:], ispin, iLsize) for j in range(Nch)] for i in range(Nch)])
    
    return np.array([idata[:,:,n] for n in range(len(idata[0,0,:]))])

def spin_proj_oct(idata, ispin, iLsize):
    idata_res = np.reshape(idata, (2,2,2*iLsize**3,2,2))
    if   (ispin == 0):
        return np.reshape(np.array([(
                        (idata_res[0,1,n,0,1] - idata_res[1,0,n,0,1]) -
                        (idata_res[0,1,n,1,0] - idata_res[1,0,n,1,0])
                        ) for n in range(2*iLsize**3)]), -1)
    elif (ispin == 1):
        return np.reshape(np.array([(
                        (idata_res[0,1,n,0,1] + idata_res[1,0,n,0,1]) +
                        (idata_res[0,1,n,1,0] + idata_res[1,0,n,1,0])
                        ) for n in range(2*iLsize**3)]), -1)
    else:
        print("\nERROR: Invalid spin.\n"); return None

def make_jk_samples (idata, Nbin):
    Nconf = len(idata[:,0])
    if (Nconf % Nbin != 0):
        print("\nERROR: Unexpected #.file, #.conf(=%d) %% #.bin(=%d) != 0.\n" % (Nconf, Nbin))
        return None
    Bsize = Nconf // Nbin
    return np.array([(np.sum(idata, axis=0)-np.sum(idata[i*Bsize:(i+1)*Bsize, :], axis=0))
                     for i in range(Nbin)]) / float(Nconf-Bsize)

def channel_projection (idata, transU):
    return np.array([np.dot(transU, np.dot(idata[n,:,:], transU.T)) for n in range(len(idata[:,0,0]))])

def calc_potential (idata):
    Lap = 1

def set_data_plot (idata, iLsize):
    Nconf     = len(idata[:,0,0,0])
    Nch       = len(idata[0,0,0,:])
    idata_res = np.reshape(idata, (Nconf, iLsize, iLsize, iLsize, 2, Nch, Nch))
    
    Lsize2    = iLsize // 2
    rdata_res = idata_res[:, :Lsize2, :Lsize2, :Lsize2, 0, :, :]
    rdata     = np.reshape(rdata_res, (Nconf, Lsize2**3, Nch, Nch))
    
    rdata_ave = np.mean(rdata, axis=0)
    if (Nconf == 1):
        rdata_err = None
    else:
        rdata_err = np.std(rdata, axis=0) * np.sqrt(Nconf-1)
    
    retR = np.reshape(np.array([[[np.sqrt(ix**2 + iy**2 + iz**2) 
                                  for ix in range(Lsize2)] for iy in range(Lsize2)] for iz in range(Lsize2)]), -1)
    
    return (retR, rdata_ave, rdata_err)

def plot_data (axe, iRcood, idata_ave, idata_err, ilabel, ilw, imarker, icolor, ixlim, iylim):
    Nch = len(idata_ave[:])
    if (idata_err is None):
        axe.errorbar(iRcood, idata_ave, fmt='.', label=ilabel,
                     lw=ilw, marker=imarker, color=icolor)
    else:
        axe.errorbar(iRcood, idata_ave, yerr=idata_err, fmt='none', label=ilabel,
                     lw=ilw, marker=imarker, ecolor=icolor)
    
    if (ixlim is not None):
        axe.set_xlim(ixlim)
    if (iylim is not None):
        axe.set_ylim(iylim)
    axe.legend(numpoints=1, loc='good')
    axe.grid(which='major',color='gray',linestyle='--')
    axe.patch.set_facecolor('white')

def plot_data_NxN (axe, iRcood, idata_ave, idata_err, ilabel, ilw, imarker, icolor):
    Nch = len(idata_ave[0,0,:])
    for i in range(Nch):
        for j in range(Nch):
            if (idata_err is None):
                plot_data(axe[i,j], iRcood, idata_ave[:,i,j], None,
                          ilabel, ilw, imarker, icolor, xlim_all, ylim_all)
            else:
                plot_data(axe[i,j], iRcood, idata_ave[:,i,j], idata_err[:,i,j],
                          ilabel, ilw, imarker, icolor, xlim_all, ylim_all)

### ============================================================ ###
### ============================================================ ###

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 2):
        exit("usage: python %s [conf list]" % basename(argv[0]))
    
    ifconf_lst = argv[1].strip()
    
    if(main(ifconf_lst) != 0):
        exit("ERROR EXIT.")
