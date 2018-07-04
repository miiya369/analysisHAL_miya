# -*- coding: utf-8 -*-

"""The module for statistics of various data."""

from numpy import mean, std, sqrt, array, sum

def make_mean_err(a_data, is_JKdata = True):
    """
    The function to calculate the mean value and error.
    
    For arguments,
    - a_data[#.data] (1-dim ndarray)
    
    return: (mean, error)
    
    Note: #.data is got from a_data.
    """
    
    if (is_JKdata):
        factor =       sqrt(len(a_data) - 1)
    else:
        factor = 1.0 / sqrt(len(a_data) - 1)
    
    return (mean(a_data), std(a_data) * factor)

def make_JKsample(a_data, a_bsize):
    """
    The function to construct the Jack-Knife samples.
    
    For arguments,
    - a_data[#.conf, #.data] (2-dim ndarray)
    
    return: data_Jack_Knife[#.bin, #.data]
    
    Note1: #.conf and #.data are got from a_data.
    Note2: #.conf % bin_size == 0 is necessary.
    """
    
    l_Nconf = len(a_data[:,0])
    if (l_Nconf % a_bsize != 0):
        print("\nERROR: Unexpected #.conf; #.conf(=%d) %% bin-size(=%d) != 0, exit.\n" %
              (l_Nconf, a_bsize)); return None
    
    l_Nbin = l_Nconf // a_bsize
    return array([(sum(a_data, axis=0)-sum(a_data[i*a_bsize:(i+1)*a_bsize, :], axis=0))
                  for i in range(l_Nbin)]) / float(l_Nconf-a_bsize)
