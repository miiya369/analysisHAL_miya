# -*- coding: utf-8 -*-

"""The module to print/output the results for the Schrodinger equation in a finite volume."""

from __future__         import print_function
from numpy              import sqrt, reshape, array
from common.statistics  import make_mean_err
from common.io_data_bin import output_bin_data_nxyz
from common.io_data_txt import output_text_data

hbar_c = 197.327053

def print_wave(a_evec, a_lat_a = None):
    """
    The function to print the wave functions.
    
    For arguments,
    - a_evec[#.conf, #.Lsite ** 3, #.state] (3-dim ndarray)
    
    Note: #.state and #.Lsite are got from a_evec.
    """
    l_Lsize = int(len(a_evec[0,:,0]) ** (1/3) + 1e-10)
    l_Nstat =     len(a_evec[0,0,:])
    
    if (l_Lsize ** 3 != len(a_evec[0,:,0])):
        print("\nERROR: Lsize is not matched, %d != %d.\n" % (l_Lsize ** 3, len(a_evec[0,:,0])))
    
    if (a_lat_a is None):
        factor = 1
    else:
        factor = a_lat_a
        
    for z in range(l_Lsize//2):
        for y in range(l_Lsize//2):
            for x in range(l_Lsize//2):
                xyz = x + l_Lsize * (y + l_Lsize * z)
                r    = sqrt(x**2 + y**2 + z**2) * factor
                tmpd = array([make_mean_err(a_evec[:,xyz,istat]) for istat in range(l_Nstat)])
                print("%lf" % r, ("%e %e "*l_Nstat) % (*tmpd.flatten(),))

def print_eigen(a_eval, a_evec, a_lat_a = None, a_eval_only = True):
    """
    The function to print the results for the Schrodinger equation in a finite volume.
    
    For arguments,
    - a_eval[#.conf,               #.state] (2-dim ndarray)
    - a_evec[#.conf, #.Lsite ** 3, #.state] (3-dim ndarray)
    
    Note: #.conf and #.state is got from a_eval.
    """
    l_Nconf = len(a_eval[:,0])
    l_Nstat = len(a_eval[0,:])
    
    if (not a_eval_only):
        print("#")
        print_wave(a_evec, a_lat_a)
    
    if (a_lat_a is None):
        factor = 1
    else:
        factor = hbar_c / a_lat_a
    
    print("#")
    for istat in range(l_Nstat):
        print("# Eigen energy[%02d] = %lf +/- %lf" % (istat, *make_mean_err(a_eval[:,istat]*factor)))
        ### For All data
        #print("# Eigen energy[%02d] =" % istat, ("%lf "*l_Nconf) % (*a_eval[:,istat]*factor,))

def output_eigen(a_ofbase, a_eval, a_evec, a_lat_a = None, a_eval_only = True):
    """
    The function to output the results for the Schrodinger equation in a finite volume.
    
    For arguments,
    - a_eval[#.conf,               #.state] (2-dim ndarray)
    - a_evec[#.conf, #.Lsite ** 3, #.state] (3-dim ndarray)
    
    Note1: #.conf and #.state are got from a_eval.
    Note2: #.Lsite is got from a_evec.
    """
    l_Nconf = len(a_eval[:,0])
    l_Nstat = len(a_eval[0,:])
    l_Lsize = int(len(a_evec[0,:,0]) ** (1/3) + 1e-10)
    
    if (l_Lsize ** 3 != len(a_evec[0,:,0])):
        print("\nERROR: Lsize is not matched, %d != %d.\n" % (l_Lsize ** 3, len(a_evec[0,:,0])))
    
    if (a_lat_a is None):
        factor = 1
    else:
        factor = hbar_c / a_lat_a
    
    if (not a_eval_only):
        for i in range(l_Nstat):
            output_bin_data_nxyz(a_ofbase+".evec.n%03d.bin" % i, 
                                 reshape(a_evec[:,:,i], (l_Nconf, l_Lsize, l_Lsize, l_Lsize)))
    
    output_text_data(a_ofbase+".eval.N%03d.txt" % l_Nstat, a_eval * factor)
