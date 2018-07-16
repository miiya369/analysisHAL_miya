# -*- coding: utf-8 -*-

"""The module for special functions."""

from numpy         import sqrt, arccos, arctan2, pi, nan, vectorize
from scipy.special import spherical_jn, spherical_yn, sph_harm

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

def sph_harm_xyz_(a_l, a_m, a_x, a_y, a_z):
    """
    The spherical harmonics in the x,y,z-coordinate.
    Defined as:
    - r     = sqrt(x^2 + y^2 + z^2)
    - theta = tan^{-1}(y/x)
    - phi   = cos^{-1}(z/r)
    
    return: sph_harm(m, n, theta, phi)
    
    Note: There are different conventions for the meanings of the 
    ....: input arguments theta and phi. In SciPy theta is the azimuthal angle 
    ....: and phi is the polar angle. It is common to see the opposite convention, 
    ....: that is, theta as the polar angle and phi as the azimuthal angle.
    """
    if (a_x == a_y == a_z == 0):
        return nan+0.0j
    l_r     = sqrt(a_x**2 + a_y**2 + a_z**2)
    l_theta = arctan2(a_y , a_x) + pi
    l_phi   = arccos (a_z / l_r)
    return sph_harm(a_m, a_l, l_theta, l_phi)

def sph_harm_xyz(a_l, a_m, a_x, a_y, a_z):
    """
    Vectorize version of sph_harm_xyz_
    """
    return vectorize(sph_harm_xyz_)(a_l, a_m, a_x, a_y, a_z)
