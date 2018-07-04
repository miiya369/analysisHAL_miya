# -*- coding: utf-8 -*-

"""The module for definition of the fit functions."""

from numpy import exp

### For curve_fit ###

def fitfunc_SW(a_x, a_p1, a_p2):
    """Square Well potential"""
    
    if (a_x < a_p1):
        return a_p2
    else:
        return 0.0

def fitfunc_Coulomb(a_x, a_p1, a_p2):
    """Coulomb potential"""
    
    # 1.43996507808 = 197.327053 / 137.035999 (= hbar c * alpha)
    
    return a_p1 * a_p2 * 1.43996507808 / a_x

def fitfunc_Coulomb_wF(a_x, a_p1, a_p2):
    """Coulomb potential w/ form factor"""
    
    # 4.27014423267 = sqrt(18.2341317678 fm^-2) (= sqrt(0.71 GeV^2))
    # 1.43996507808 = 197.327053 / 137.035999   (= hbar c * alpha)
    
    sqrt_b_r = 4.27014423267 * a_x
    factor   = (2.0 - exp(-sqrt_b_r) * (2.0 + sqrt_b_r)) / 2.0
    
    return a_p1 * a_p2 * 1.43996507808 / a_x * factor

def fitfunc_1G(a_x, a_p1, a_p2):
    """one Gaussian"""
    
    return (a_p1 * exp(-pow(a_x/a_p2, 2)))

def fitfunc_2G(a_x, a_p1, a_p2, a_p3, a_p4):
    """2-ranges Gaussian"""
    
    return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)))

def fitfunc_3G(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6):
    """3-ranges Gaussian"""
    
    return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
            a_p5 * exp(-pow(a_x/a_p6, 2)))

def fitfunc_4G(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8):
    """4-ranges Gaussian"""
    
    return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
            a_p5 * exp(-pow(a_x/a_p6, 2)) + a_p7 * exp(-pow(a_x/a_p8, 2)))

def fitfunc_12G(a_x, 
                a_p1 , a_p2 , a_p3 , a_p4 , a_p5 , a_p6 , a_p7 , a_p8,
                a_p9 , a_p10, a_p11, a_p12, a_p13, a_p14, a_p15, a_p16,
                a_p17, a_p18, a_p19, a_p20, a_p21, a_p22, a_p23, a_p24):
    """12-ranges Gaussian"""
    
    return (a_p1  * exp(-pow(a_x/a_p2 , 2)) + a_p3  * exp(-pow(a_x/a_p4 , 2)) +
            a_p5  * exp(-pow(a_x/a_p6 , 2)) + a_p7  * exp(-pow(a_x/a_p8 , 2)) +
            a_p9  * exp(-pow(a_x/a_p10, 2)) + a_p11 * exp(-pow(a_x/a_p12, 2)) +
            a_p13 * exp(-pow(a_x/a_p14, 2)) + a_p15 * exp(-pow(a_x/a_p16, 2)) +
            a_p17 * exp(-pow(a_x/a_p18, 2)) + a_p19 * exp(-pow(a_x/a_p20, 2)) +
            a_p21 * exp(-pow(a_x/a_p22, 2)) + a_p23 * exp(-pow(a_x/a_p24, 2)))

def fitfunc_1SG(a_x, a_p1, a_p2, a_p3):
    """one shifted-Gaussian"""
    
    return (a_p1 * exp(-pow((a_x-a_p2)/a_p3, 2)))

def fitfunc_2SG(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6):
    """2-ranges shifted-Gaussian"""
    
    return (a_p1 * exp(-pow((a_x-a_p2)/a_p3, 2)) + a_p4 * exp(-pow((a_x-a_p5)/a_p6, 2)))

def fitfunc_1Exp(a_x, a_p1, a_p2):
    """one Exponential"""
    
    return (a_p1 * exp(-a_p2 * a_x))

def fitfunc_1Y(a_x, a_p1, a_p2, a_p3):
    """one Yukawa"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * (1.0 - exp(-a_p2 * a_x**2)) * exp(-a_p3 * a_x) / a_x)
    else:
        return 0.0
    
def fitfunc_2Y(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6):
    """2-ranges Yukawa"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * (1.0 - exp(-a_p2 * a_x**2)) * exp(-a_p3 * a_x) / a_x +
                a_p4 * (1.0 - exp(-a_p5 * a_x**2)) * exp(-a_p6 * a_x) / a_x)
    else:
        return 0.0

def fitfunc_1Ysq(a_x, a_p1, a_p2, a_p3):
    """one Yukawa square"""

#    if (a_x != 0.0):
    if (True):
        return (a_p1 * (1.0 - exp(-a_p2 * a_x**2))**2 * exp(-2.0 * a_p3 * a_x) / a_x**2)
    else:
        return 0.0

def fitfunc_2Ysq(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6):
    """2-ranges Yukawa square"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * (1.0 - exp(-a_p2 * a_x**2))**2 * exp(-2.0 * a_p3 * a_x) / a_x**2 +
                a_p4 * (1.0 - exp(-a_p5 * a_x**2))**2 * exp(-2.0 * a_p6 * a_x) / a_x**2)
    else:
        return 0.0

def fitfunc_1Ytns(a_x, a_p1, a_p2, a_p3):
    """one Yukawa tensor"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * (1.0 - exp(-a_p2 * a_x**2))**2 * (1.0 + 3.0/(a_p3*a_x) + 3.0/(a_p3*a_x)**2)
                * exp(-a_p3 * a_x) / a_x)
    else:
        return 0.0

def fitfunc_2Ytns(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6):
    """2-ranges Yukawa tensor"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * (1.0 - exp(-a_p2 * a_x**2))**2 * (1.0 + 3.0/(a_p3*a_x) + 3.0/(a_p3*a_x)**2)
                * exp(-a_p3 * a_x) / a_x +
                a_p4 * (1.0 - exp(-a_p5 * a_x**2))**2 * (1.0 + 3.0/(a_p6*a_x) + 3.0/(a_p6*a_x)**2)
                * exp(-a_p6 * a_x) / a_x)
    else:
        return 0.0

def fitfunc_1G1Y(a_x, a_p1, a_p2, a_p3, a_p4, a_p5):
    """1-ranges Gaussian + one Yukawa"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * (1.0 - exp(-a_p4 * a_x**2)) * exp(-a_p5  * a_x) / a_x)
    else:
        return  a_p1

def fitfunc_2G1Y(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7):
    """2-ranges Gaussian + one Yukawa"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
                a_p5 * (1.0 - exp(-a_p6 * a_x**2)) * exp(-a_p7  * a_x) / a_x)
    else:
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)))

def fitfunc_2G2Y(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9, a_p10):
    """2-ranges Gaussian + 2-ranges Yukawa"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
                a_p5 * (1.0 - exp(-a_p6 * a_x**2)) * exp(-a_p7  * a_x) / a_x +
                a_p8 * (1.0 - exp(-a_p9 * a_x**2)) * exp(-a_p10 * a_x) / a_x)
    else:
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)))

def fitfunc_2G1Ysq(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7):
    """2-ranges Gaussian + one Yukawa square"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
                a_p5 * (1.0 - exp(-a_p6 * a_x**2))**2 * exp(-2.0 * a_p7 * a_x) / a_x**2)
    else:
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)))

def fitfunc_3G1Ysq(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9):
    """3-ranges Gaussian + one Yukawa square"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
                a_p5 * exp(-pow(a_x/a_p6, 2)) +
                a_p7 * (1.0 - exp(-a_p8 * a_x**2))**2 * exp(-2.0 * a_p9 * a_x) / a_x**2)
    else:
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
                a_p5 * exp(-pow(a_x/a_p6, 2)))

def fitfunc_2G1Y1Ysq(a_x, a_p1, a_p2, a_p3, a_p4, a_p5, a_p6, a_p7, a_p8, a_p9, a_p10):
    """2-ranges Gaussian + one Yukawa + one Yukawa square"""
    
#    if (a_x != 0.0):
    if (True):
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)) +
                a_p5 * (1.0 - exp(-a_p6 * a_x**2))    * exp(-      a_p7  * a_x) / a_x +
                a_p8 * (1.0 - exp(-a_p9 * a_x**2))**2 * exp(-2.0 * a_p10 * a_x) / a_x**2)
    else:
        return (a_p1 * exp(-pow(a_x/a_p2, 2)) + a_p3 * exp(-pow(a_x/a_p4, 2)))
