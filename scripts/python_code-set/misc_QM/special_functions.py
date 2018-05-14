# -*- coding: utf-8 -*-

"""The module for special functions."""

from scipy.special import spherical_jn, spherical_yn, riccati_jn, riccati_yn
from numpy         import vectorize

def spherical_hn1(a_n, a_z, a_d=False):
    """
    The spherical Hankel function of the first kind.
    
    return: jn + i yn
    """
    return spherical_jn(a_n, a_z, a_d) + 1.0j * spherical_yn(a_n, a_z, a_d)

def spherical_hn2(a_n, a_z, a_d=False):
    """
    The spherical Hankel function of the second kind.
    
    return: jn - i yn
    """
    return spherical_jn(a_n, a_z, a_d) - 1.0j * spherical_yn(a_n, a_z, a_d)

def riccati_jn_wrap(a_n, a_z, a_d=False):
    """
    The wrapper function for the Riccati-Bessel function.
    
    return: z * jn(z) (or its derivative)
    """
    if (a_d):
        return riccati_jn(a_n, a_z)[1][a_n]
    else:
        return riccati_jn(a_n, a_z)[0][a_n]

def riccati_yn_wrap(a_n, a_z, a_d=False):
    """
    The wrapper function for the Riccati-Neumann function.
    
    return: z * yn(z) (or its derivative)
    """
    if (a_d):
        return riccati_yn(a_n, a_z)[1][a_n]
    else:
        return riccati_yn(a_n, a_z)[0][a_n]

def riccati_jn_vec(a_n, a_z, a_d=False):
    """
    The vectorize version for the Riccati-Bessel function.
    
    return: z * jn(z) (or its derivative)
    """
    return vectorize(riccati_jn_wrap)(a_n, a_z, a_d)

def riccati_yn_vec(a_n, a_z, a_d=False):
    """
    The vectorize version for the Riccati-Neumann function.
    
    return: z * yn(z) (or its derivative)
    """
    return vectorize(riccati_yn_wrap)(a_n, a_z, a_d)

def riccati_hn1(a_n, a_z, a_d=False):
    """
    The Riccati-Hankel function of the first kind.
    Defined as z * spherical_h1(z).
    
    return: jn + i yn
    """
    return vectorize(riccati_jn_wrap)(a_n, a_z, a_d) + 1.0j * vectorize(riccati_yn_wrap)(a_n, a_z, a_d)

def riccati_hn2(a_n, a_z, a_d=False):
    """
    The Riccati-Hankel function of the second kind.
    Defined as z * spherical_h2(z).
    
    return: jn - i yn
    """
    return vectorize(riccati_jn_wrap)(a_n, a_z, a_d) - 1.0j * vectorize(riccati_yn_wrap)(a_n, a_z, a_d)
