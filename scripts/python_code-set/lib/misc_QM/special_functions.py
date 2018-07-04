# -*- coding: utf-8 -*-

"""The module for special functions."""

from scipy.special import spherical_jn, spherical_yn

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

def riccati_jn(a_n, a_z, a_d=False):
    """
    The Riccati-Bessel function.
    
    return: z * jn(z) (or its derivative)
    """
    if (a_d):
        return spherical_jn(a_n, a_z) + a_z * spherical_jn(a_n, a_z, True)
    else:
        return a_z * spherical_jn(a_n, a_z)

def riccati_yn(a_n, a_z, a_d=False):
    """
    The Riccati-Neumann function.
    
    return: z * yn(z) (or its derivative)
    """
    if (a_d):
        return spherical_yn(a_n, a_z) + a_z * spherical_yn(a_n, a_z, True)
    else:
        return a_z * spherical_yn(a_n, a_z)

def riccati_hn1(a_n, a_z, a_d=False):
    """
    The Riccati-Hankel function of the first kind.
    Defined as z * spherical_h1(z).
    
    return: jn + i yn
    """
    return riccati_jn(a_n, a_z, a_d) + 1.0j * riccati_yn(a_n, a_z, a_d)

def riccati_hn2(a_n, a_z, a_d=False):
    """
    The Riccati-Hankel function of the second kind.
    Defined as z * spherical_h2(z).
    
    return: jn - i yn
    """
    return riccati_jn(a_n, a_z, a_d) - 1.0j * riccati_yn(a_n, a_z, a_d)
