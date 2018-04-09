# -*- coding: utf-8 -*-

"""The module to print/output the results for the Schrodinger equation in a finite volume."""

from __future__         import print_function
from numpy              import sqrt, reshape, array
from common.statistics  import make_mean_err
from common.io_data_bin import output_bin_data_nxyz
from common.io_data_txt import output_text_data

hbar_c = 197.327053

def print_wave(a_EigVec, a_lat_a = None):
    """
    The function to print the wave functions.
    
    For arguments,
    - a_EigVec[#.conf, #.Lsite ** 3, #.state] (3-dim ndarray)
    
    Note: #.state and #.Lsite are got from a_EigVec.
    """
    
    l_Lsize  = int(len(a_EigVec[0,:,0]) ** (1/3) + 1e-10)
    l_Nstate =     len(a_EigVec[0,0,:])
    
    if (l_Lsize ** 3 != len(a_EigVec[0,:,0])):
        print("\nERROR: Lsize is not matched, %d != %d.\n" % (l_Lsize ** 3, len(a_EigVec[0,:,0])))
    
    if (a_lat_a is None):
        factor = 1
    else:
        factor = a_lat_a
        
    for z in range(l_Lsize//2):
        for y in range(l_Lsize//2):
            for x in range(l_Lsize//2):
                ixyz = x + l_Lsize * (y + l_Lsize * z)
                ir   = sqrt(x**2 + y**2 + z**2) * factor
                print("%f" % ir, end="")
                for istate in range(l_Nstate):
                    mean, err = make_mean_err(a_EigVec[:,ixyz,istate])
                    print(" %e %e" % (mean, err), end="")
                print()

def print_results(a_EigVal, a_EigVec, a_lat_a = None, a_EigVal_Only = True):
    """
    The function to print the results for the Schrodinger equation in a finite volume.
    
    For arguments,
    - a_EigVal[#.conf,               #.state] (2-dim ndarray)
    - a_EigVec[#.conf, #.Lsite ** 3, #.state] (3-dim ndarray)
    
    Note: #.state is got from a_EigVal.
    """
    
    l_Nstate = len(a_EigVal[0,:])
    
    if (a_lat_a is None):
        factor = 1
    else:
        factor = hbar_c / a_lat_a
    
    if (not a_EigVal_Only):
        print("#")
        print_wave(a_EigVec, a_lat_a)
    
    print("#")
    for istate in range(l_Nstate):
        mean, err = make_mean_err(a_EigVal[:,istate] * factor)
        print("# Eigen energy[%02d] = %lf +/- %lf" % 
              (istate, mean, err))

def output_results(a_ofbase, a_EigVal, a_EigVec, a_lat_a = None, a_EigVal_Only = True):
    """
    The function to output the results for the Schrodinger equation in a finite volume.
    
    For arguments,
    - a_EigVal[#.conf,               #.state] (2-dim ndarray)
    - a_EigVec[#.conf, #.Lsite ** 3, #.state] (3-dim ndarray)
    
    Note1: #.conf and #.state are got from a_EigVal.
    Note2: #.Lsite is got from a_EigVec.
    """
    
    l_Nconf  = len(a_EigVal[:,0])
    l_Nstate = len(a_EigVal[0,:])
    l_Lsize  = int(len(a_EigVec[0,:,0]) ** (1/3) + 1e-10)
    
    if (l_Lsize ** 3 != len(a_EigVec[0,:,0])):
        print("\nERROR: Lsize is not matched, %d != %d.\n" % (l_Lsize ** 3, len(a_EigVec[0,:,0])))
    
    if (a_lat_a is None):
        factor = 1
    else:
        factor = hbar_c / a_lat_a
    
    if (not a_EigVal_Only):
        for i in range(l_Nstate):
            output_bin_data_nxyz(a_ofbase+".Evec.n%03d.bin" % i, 
                                 reshape(a_EigVec[:,:,i], (l_Nconf, l_Lsize, l_Lsize, l_Lsize)))
    
    output_text_data(a_ofbase+".Eval.N%03d.txt" % l_Nstate, a_EigVal * factor)
