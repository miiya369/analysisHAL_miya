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

def sph_harm_A1_xyz(a_l, a_m, a_x, a_y, a_z):
    """
    A1 projected spherical harmonics (a_l < 9)
    """
    if   (a_l == 0):
        return 1.0 / (2.0*sqrt(pi))
    
    elif (a_l == 1):
        return a_x * a_y * a_z * 0.0
    elif (a_l == 2):
        return a_x * a_y * a_z * 0.0
    elif (a_l == 3):
        return a_x * a_y * a_z * 0.0
    
    elif (a_l == 4):
        X2 = a_x**2; Y2 = a_y**2; Z2 = a_z**2
        X4 = a_x**4; Y4 = a_y**4; Z4 = a_z**4
        nume =  X4 + Y4 + Z4 - 3*(X2*Y2 + Y2*Z2 + Z2*X2)
        deno = (X2 + Y2 + Z2)**2
        if   (    a_m  == 0):
            return (7.0/(8.0*sqrt(pi)))      * (float(nume)/float(deno))
        elif (abs(a_m) == 4):
            return (sqrt(35.0/(2.0*pi))/8.0) * (float(nume)/float(deno))
        else:
            return a_x * a_y * a_z * 0.0
    
    elif (a_l == 5):
        return a_x * a_y * a_z * 0.0
    
    elif (a_l == 6):
        X2 = a_x**2; Y2 = a_y**2; Z2 = a_z**2
        X4 = a_x**4; Y4 = a_y**4; Z4 = a_z**4
        X6 = a_x**6; Y6 = a_y**6; Z6 = a_z**6
        nume = 2*(X6 + Y6 + Z6) - 15*(X2*Y4 + Y2*Z4 + Z2*X4 +
                                      X4*Y2 + Y4*Z2 + Z4*X2) + 180*X2*Y2*Z2
        deno =   (X2 + Y2 + Z2)**3
        if   (    a_m  == 0):
            return  (sqrt(13.0/     pi )/32.0) * (float(nume)/float(deno))
        elif (abs(a_m) == 4):
            return -(sqrt(91.0/(2.0*pi))/32.0) * (float(nume)/float(deno))
        else:
            return a_x * a_y * a_z * 0.0
    
    elif (a_l == 7):
        return a_x * a_y * a_z * 0.0
    
    elif (a_l == 8):
        X2 = a_x**2; Y2 = a_y**2; Z2 = a_z**2
        X4 = a_x**4; Y4 = a_y**4; Z4 = a_z**4
        X6 = a_x**6; Y6 = a_y**6; Z6 = a_z**6
        X8 = a_x**8; Y8 = a_y**8; Z8 = a_z**8
        nume =  X8 + Y8 + Z8 - 14*(X2*Y6 + Y2*Z6 + Z2*X6 +
                                   X6*Y2 + Y6*Z2 + Z6*X2) + 35*(X4*Y4 + Y4*Z4 + Z4*X4)
        deno = (X2 + Y2 + Z2)**4
        if   (    a_m  == 0):
            return (sqrt(   17.0/pi) *33.0/128.0) * (float(nume)/float(deno))
        elif (abs(a_m) == 4):
            return (sqrt( 1309.0/(2.0*pi))/ 64.0) * (float(nume)/float(deno))
        elif (abs(a_m) == 8):
            return (sqrt(12155.0/(2.0*pi))/128.0) * (float(nume)/float(deno))
        else:
            return a_x * a_y * a_z * 0.0        
    
    else:
        print("\nERROR: Invalid (or have not been implemented) " +
              "angular momentum L=%d.\n" % a_l)
        return None
