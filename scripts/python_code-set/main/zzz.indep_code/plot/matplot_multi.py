#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc, patches
rc('text' , usetex=True)
rc('font' ,**{'family'    : 'Times New Roman',
              'weight'    : 'bold',
              'size'      : 12})
rc('xtick',**{'labelsize' : 14})
rc('ytick',**{'labelsize' : 14})
rc('axes' ,**{'labelsize' : 14})
rc('axes' ,**{'linewidth' : 2})

### =========================== Main =========================== ###

def main (ifpath, ofname):
    if (ofname == 'None'):
        plt.ion()
    
    plot_sub2(ifpath)
    
    if (ofname == 'None'):
        plt.show(); print("# push Enter to end"); input('# ') #raw_input('# ') # For v.2
    else:
        print("Output file: %s" % ofname)
        plt.savefig(ofname)
    
    return 0

### ============================================================ ###

def plot_sub3(ifpath): # For PoS HADRON 2017 (Phase shift)
    ifname_base = "%s/ens%d.gnu"
    ylabel      = ('[degrees]','[degrees]','')
    xlabel      = (r'$E_{\mathrm{CM}} - m_{\Lambda_cN}$ [MeV]', 
                   r'$E_{\mathrm{CM}} - m_{\Sigma_cN}$ [MeV]', 
                   r'$E_{\mathrm{CM}} - m_{\Sigma_cN}$ [MeV]')
    labels      = (r'$m_\pi \simeq 700$ MeV', r'$m_\pi \simeq 570$ MeV', r'$m_\pi \simeq 410$ MeV')
    titles      = (r'(a) $^3S_1$ $\Lambda_c N$ phase shift $\bar{\delta}_{\Lambda_cN}$',
                   r'(b) $^3S_1$ $\Sigma_c N$ phase shift $\bar{\delta}_{\Sigma_cN}$',
                   r'(c) Inelasticity of the scattering $|S_{\Lambda_cN}^{~\Lambda_cN}|$')
    FigSizeBase = [15, 4]
    FigSize     = [[0.05, 0.15, 0.25, 0.75], [0.35, 0.15, 0.25, 0.75], 
                   [0.66, 0.15, 0.25, 0.75], [0.72, 0.45, 0.18, 0.30]]
    yrange      = [[-20, 50], [-20, 50], [0.9, 1.01], [0.99, 1.00]]
    xrange      = [[0, 180], [0, 80], [0, 80], [0, 40]]
    yaxs_cord   = (-0.1,0.5)
    xaxs_cord   = (0.5,-0.12)
    Legend_pls  = [1.0, 1.0]
    colors      = ['blue','lime','red']
    
    th = (95.7374,119.135,140.806)
    
    fig = plt.figure(figsize=FigSizeBase)
    fig.patch.set_facecolor('white')
    
    axe = [None for i in range(4)]
    for i in range(4):
        axe[i] = fig.add_axes(FigSize[i])
        axe[i].patch.set_facecolor('white')
        if (i != 3):
            axe[i].set_xlabel(xlabel[i])
            axe[i].xaxis.set_label_coords(xaxs_cord[0],xaxs_cord[1])
            axe[i].set_title(titles[i])
            axe[i].set_ylabel(ylabel[i])
            axe[i].yaxis.set_label_coords(yaxs_cord[0],yaxs_cord[1])
    
    for ens in (1,2,3):
        ifname = ifname_base % (ifpath, ens)
        print("Load: %s" % ifname); idata = np.loadtxt(ifname)
        axe[0].errorbar(idata[ens-1::3,0], idata[ens-1::3,1], yerr=idata[ens-1::3,2], fmt='none', #fmt=None, 
                        marker=None, label=labels[ens-1], lw=2, ecolor=colors[ens-1], capsize=0)
        axe[1].errorbar(idata[ens-1::2,0]-th[ens-1], idata[ens-1::2,3], 
                        yerr=idata[ens-1::2,4], fmt='none', #fmt=None, 
                        marker=None, label=labels[ens-1], lw=3, ecolor=colors[ens-1], capsize=0)
        axe[2].errorbar(idata[ens-1::2,0]-th[ens-1], idata[ens-1::2,7], 
                        yerr=idata[ens-1::2,8], fmt='none', #fmt=None, 
                        marker=None, label=labels[ens-1], lw=3, ecolor=colors[ens-1], capsize=0)
        axe[3].errorbar(idata[ens-1::2,0]-th[ens-1], idata[ens-1::2,7], 
                        yerr=idata[ens-1::2,8], fmt='none', #fmt=None, 
                        marker=None, label=None, lw=3, ecolor=colors[ens-1], capsize=0)
    
    for i in range(4):
        axe[i].set_ylim(yrange[i]); axe[i].set_xlim(xrange[i])
        axe[i].axhline(0.0, lw=2, color='black')
        axe[i].grid(which='major',color='gray',linestyle='--')
        if (i != 3):
            axe[i].legend(numpoints=1, loc='upper right', bbox_to_anchor=Legend_pls)
    axe[1].legend(numpoints=1, loc='lower right')#, bbox_to_anchor=Legend_pls)
    axe[2].legend(numpoints=1, loc='lower right')#, bbox_to_anchor=Legend_pls)
    axe[2].axhline(1.0, lw=2, color='black')

def plot_sub2(ifpath): # For PoS HADRON 2017 (Potential)
    ifname_base = "%s/L%sN_spin%d_ens%d_t%02d_CCP%d%d_err"
    ylabel      = r'$V_C^{~C^\prime} (r)$ [MeV]'
    xlabel      = r'$r$ [fm]'
    quark       = "c"
    spin        = 1
    itime       = (13, 11, 9)
    plot_ch     = ((0,0), (1,1), (0,1), (1,0))
    labels      = (r'$m_\pi \simeq 700$ MeV', r'$m_\pi \simeq 570$ MeV', r'$m_\pi \simeq 410$ MeV')
    titles      = (r'(a) $V_{\Lambda_cN}^{~\Lambda_cN}$',
                   r'(b) $V_{\Sigma_cN}^{~\Sigma_cN}$',
                   r'(c) Average of the off-diagonal elements')
    FigSizeBase = [15, 4]
    FigSize     = [[[0.07, 0.15, 0.25, 0.75], [0.13, 0.35, 0.18, 0.5]],
                   [[0.37, 0.15, 0.25, 0.75], [0.43, 0.35, 0.18, 0.5]],
                   [[0.67, 0.15, 0.25, 0.75], [0.73, 0.35, 0.18, 0.5]]]
    yrange      = [[[-200, 2000], [-50, 60]],
                   [[-200, 2000], [-50, 60]],
                   [[-200, 2000], [-50, 60]]]
    xrange      = [[-0.1, 2.0], [0.0, 2.0]]
    yaxs_cord   = (-0.15,0.5)
    xaxs_cord   = (0.5,-0.12)
    Legend_pls  = [1.0, 1.0]
    colors      = ['blue','lime','red']
    
    fig = plt.figure(figsize=FigSizeBase)
    fig.patch.set_facecolor('white')
    
    axe = [[None for j in range(2)] for i in range(3)]
    for i in range(3):
        for j in range(2):
            axe[i][j] = fig.add_axes(FigSize[i][j])
            axe[i][j].patch.set_facecolor('white')
        axe[i][0].set_xlabel(xlabel)
        axe[i][0].xaxis.set_label_coords(xaxs_cord[0],xaxs_cord[1])
        axe[i][0].set_title(titles[i])
    axe[0][0].set_ylabel(ylabel)
    axe[0][0].yaxis.set_label_coords(yaxs_cord[0],yaxs_cord[1])
    
    for ens in (1,2,3):
        for i in range(3):
            ifname = ifname_base % (ifpath, quark, spin, ens, itime[ens-1], plot_ch[i][0], plot_ch[i][1])
            print("Load: %s" % ifname); idata = np.loadtxt(ifname)
            if (i == 2):
                ifname = ifname_base % (ifpath, quark, spin, ens, itime[ens-1], plot_ch[i+1][0], plot_ch[i+1][1])
                print("Load: %s" % ifname); idata += np.loadtxt(ifname); idata /= 2.0
            
            axe[i][0].errorbar(idata[:,0], idata[:,1], yerr=idata[:,2], fmt='none', #fmt=None, 
                               marker=None, label=None         , lw=4, ecolor=colors[ens-1], capthick=4, capsize=2)
            axe[i][1].errorbar(idata[:,0], idata[:,1], yerr=idata[:,2], fmt='none', #fmt=None, 
                               marker=None, label=labels[ens-1], lw=2, ecolor=colors[ens-1],             capsize=0)
    
    for i in range(3):
        for j in range(2):
            axe[i][j].set_ylim(yrange[i][j]); axe[i][j].set_xlim(xrange[j])
            axe[i][j].axhline(0.0, lw=2, color='black')
            axe[i][j].grid(which='major',color='gray',linestyle='--')
        #axe[i][0].set_yticks([i for i in range(yrange[i][0][0], yrange[i][0][1]+1,200)])
        #axe[i][1].set_yticks([i for i in range(yrange[i][1][0], yrange[i][1][1]+1, 20)])
        axe[i][1].legend(numpoints=1, loc='upper right', bbox_to_anchor=Legend_pls)

def plot_sub1(ifpath): # For 2x2 potential 
    Nch        = 2
    ylabel     = r'$^3S_1$ effective-central potential $V(r)$ [MeV]'
    xlabel     = r'$r$ [fm]'
    titles     = (('27 (snk) - 27 (src)', '27 (src) - 8s (src)'),
                  ('8s (snk) - 27 (src)', '8s (snk) - 8s (src)'))
    FigSize    = [14, 10]
    fig_space  = (0.2, 0.2)
    yaxs_cord  = (-0.25,-0.1)
    xaxs_cord  = (0.5,-0.13)
    xlim_all   = None
    ylim_all   = None
    Legend_pls = [1.0, 1.0]
    colors     = ['cyan','magenta','blue','lime','red','yellow','black']
    
    fig, axe = plt.subplots(nrows=Nch, ncols=Nch, figsize=FigSize)
    plt.subplots_adjust(wspace=fig_space[0], hspace=fig_space[1])
    fig.patch.set_facecolor('white')
    axe[0,0].set_ylabel(ylabel)
    axe[0,0].yaxis.set_label_coords(yaxs_cord[0], yaxs_cord[1])
    
    for i in range(Nch):
        for j in range(Nch):
            print("Load: %s" % ifpath); idata = np.loadtxt(ifpath)
            axe[i,j].patch.set_facecolor('white')
            axe[i,j].set_title(r'%s' % titles[i][j])
            
            axe[i,j].errorbar(idata[:,0], idata[:,1], yerr=idata[:,2], fmt='none', #fmt=None, 
                              marker=None, label=None, lw=3, ecolor='red', capthick=3, capsize=3)
            
            if (ylim_all is not None):
                axe[i,j].set_ylim(ylim_all)
            if (xlim_all is not None):
                axe[i,j].set_xlim(xlim_all)
            
            axe[i,j].axhline(0.0, lw=2, color='black')
            axe[i,j].grid(which='major',color='gray',linestyle='--')
            axe[i,j].legend(numpoints=1, loc='upper right', bbox_to_anchor=Legend_pls)
            axe[i,j].set_xlabel(xlabel)
            axe[i,j].xaxis.set_label_coords(xaxs_cord[0], xaxs_cord[1])
### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 3):
        exit("usage: python %s [idir] [ofname]" % basename(argv[0]))
    
    if (main(argv[1], argv[2]) != 0):
        exit("ERROR EXIT.")
