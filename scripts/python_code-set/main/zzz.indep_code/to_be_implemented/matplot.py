#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Sun May 28 00:27:38 JST 2017
### Brief : Plot the data using matplotlib

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc
rc('text' , usetex=True)
rc('font' ,**{'family'    :'sans-serif',
              'weight'    : 'bold',
              'size'      : 14,
              'sans-serif':['Helvetica']})
rc('xtick',**{'labelsize' : 18})
rc('ytick',**{'labelsize' : 18})
rc('axes' ,**{'labelsize' : 18})
rc('axes' ,**{'linewidth' : 2})

Fig1_size = [0.15, 0.15, 0.8, 0.8]
Fig2_size = [0.35, 0.35, 0.57, 0.57]

Legend_pls = [1.0, 1.0]

yrange1 = [-500, 3500]
yrange2 = [-30, 40]
xrange  = [0, 2.0]

ylabel = r'$^1S_0$ $\Lambda_c N$ central potential $V_C(r)$ [MeV]'
xlabel = r'$r$ [fm]'

labels = [r'$m_\pi \sim 700$ MeV',
          r'$m_\pi \sim 570$ MeV',
          r'$m_\pi \sim 410$ MeV']

plot_idx = [0,1,2]

ofname = None
#ofname = "tmp.eps"

def main(ARGC, ARGV):
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    
    ax1 = fig.add_axes(Fig1_size)
    ax1.patch.set_facecolor('white')
    ax2 = fig.add_axes(Fig2_size, sharex=ax1)
    ax2.patch.set_facecolor('white')
    
    for i in range(1, ARGC):
        print("Load: %s" % ARGV[i])
        idata = np.loadtxt(ARGV[i].strip())
        ax1.errorbar(idata[:,plot_idx[0]], idata[:,plot_idx[1]], yerr=idata[:,plot_idx[2]],
                     fmt='.', label=labels[i-1])
        ax2.errorbar(idata[:,plot_idx[0]], idata[:,plot_idx[1]], yerr=idata[:,plot_idx[2]],
                     fmt='.', label=labels[i-1])
    
    ax1.set_ylim(yrange1)
    ax1.set_xlim(xrange )
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax1.axhline(0.0, lw=2, color='black')
    ax1.grid(which='major',color='gray',linestyle='--')
    
    ax2.set_ylim(yrange2)
    ax2.axhline(0.0, lw=2, color='black')
    ax2.legend(numpoints=1, loc='upper right', bbox_to_anchor=Legend_pls)
    ax2.grid(which='major',color='gray',linestyle='--')

if (__name__ == "__main__"):
    from sys     import argv; argc = len(argv)
    from os.path import basename
    
    if (argc < 2):
        print("usage: python %s [text data(s)]" % basename(argv[0])); quit()
    
    plt.ion (); main(argc, argv)
    
    if (ofname is None):
        plt.show(); print("@@@ push Enter to end"),; raw_input(': ')
    else:
        print("Output file: %s" % ofname)
        plt.savefig(ofname)
