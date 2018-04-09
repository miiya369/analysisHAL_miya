# -*- coding: utf-8 -*-

"""The module for statistics of various data."""

from numpy import mean, std, sqrt, array, sum

def make_mean_err(a_Data, is_JKdata = True):
    """
    The function to calculate the mean value and error.
    
    For arguments,
    - a_Data[#.data] (1-dim ndarray)
    
    return: (mean, error)
    
    Note: #.data is got from a_Data.
    """
    
    if (is_JKdata):
        factor =       sqrt(len(a_Data) - 1)
    else:
        factor = 1.0 / sqrt(len(a_Data) - 1)
    
    return (mean(a_Data), std(a_Data) * factor)

def make_JKsample(a_Data, a_Bsize):
    """
    The function to construct the Jack-Knife samples.
    
    For arguments,
    - a_Data[#.conf, #.data] (2-dim ndarray)
    
    return: Data_Jack_Knife[#.bin, #.data]
    
    Note1: #.conf and #.data are got from a_Data.
    Note2: #.conf % bin_size == 0 is necessary.
    """
    
    l_Nconf = len(a_Data[:,0])
    if (l_Nconf % a_Bsize != 0):
        print("\nERROR: Unexpected #.conf; #.conf(=%d) %% bin-size(=%d) != 0, exit.\n" %
              (l_Nconf, a_Bsize)); return None
    
    l_Nbin = l_Nconf // a_Bsize
    return array([(sum(a_Data, axis=0)-sum(a_Data[i*a_Bsize:(i+1)*a_Bsize, :], axis=0))
                  for i in range(l_Nbin)]) / float(l_Nconf-a_Bsize)
