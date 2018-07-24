# -*- coding: utf-8 -*-

"""The module to calculate phase shifts by Luscher's finite volume formula."""

from math            import factorial, ceil
from numpy           import pi, exp, sqrt
from scipy.integrate import quad

def sum_0_inf(a_f, precision = 10e-9):
    """
    Sub function to calculate sum_{x=0}^inf f(x)
    
    return: sum_{x=0}^inf f(x)
    """    
    l_MaxIter   = 1000000
    l_IterBlock = 100
    
    ret         = 0.0
    is_converge = False
    for i in range(l_MaxIter//l_IterBlock):
        for x in range(l_IterBlock):
            ret += a_f(x + l_IterBlock * i)
        if (abs(a_f(l_IterBlock * (i+1))) < precision):
            ret += a_f(l_IterBlock * (i+1))
            is_converge = True
            break
    
    if (not is_converge):
        print("ERROR in sum_0_inf: Not converge.")
        ret = None
    
    return ret

def sum_1_Z3_old(a_f, precision = 10e-9):
    """
    Sub function to calculate sum_{n in Z^3 (without n=0)} f(|n|^2)
    
    return: sum_{n in Z^3 (without n=0)} f(|n|^2)
    """
    l_Max_n = 10000
    
    ret         = 0.0
    is_converge = False
    for n_in in range(1, l_Max_n+1):
        for n1 in range(-n_in, n_in+1):
            for n2 in range(-n_in, n_in+1):
                for n3 in range(-n_in, n_in+1):
                    if (n1**2 + n2**2 + n3**2 != n_in):
                        continue
                    ret += a_f(n_in)
        
        if (abs(a_f(n_in)) < precision):
            is_converge = True
            break
    
    if (not is_converge):
        print("ERROR in sum_1_Z3: Not converge.")
        ret = None
    
    return ret

def sum_1_Z3(a_f, precision = 10e-9):
    """
    Sub function to calculate sum_{n in Z^3 (without n=0)} f(|n|^2)
    
    return: sum_{n in Z^3 (without n=0)} f(|n|^2)
    """
    l_Max_n = 10000
    
    ret         = 0.0
    is_converge = False
    for n_in in range(1, l_Max_n+1):
        ret += a_f(n_in) * num_point_r2_rev2(n_in)
        
        if (abs(a_f(n_in)) < precision):
            is_converge = True
            break
    
    if (not is_converge):
        print("ERROR in sum_1_Z3: Not converge.")
        ret = None
    
    return ret

def num_point_r2(a_r2):
    """
    Return the number of r^2 which satisfies x^2 + y^2 + z^2 = r^2, where x, y, z are in Z.
    """
    ret = 0
    for x in range(-a_r2, a_r2+1):
        for y in range(-a_r2, a_r2+1):
            for z in range(-a_r2, a_r2+1):
                if (x**2 + y**2 + z**2 == a_r2):
                    ret += 1
    return ret

def num_point_r2_rev(a_r2):
    """
    Return the number of r^2 which satisfies x^2 + y^2 + z^2 = r^2, where x, y, z are in Z.
    ~~~ Modified version ~~~
    """
    ret = 0
    if (a_r2 == 0):
        return 1
    
    if(int(sqrt(a_r2/1.0))**2 == a_r2/1.0):
        ret += 6
    if(int(sqrt(a_r2/2.0))**2 == a_r2/2.0):
        ret += 12
    if(int(sqrt(a_r2/3.0))**2 == a_r2/3.0):
        ret += 8
    
    tmp_i = 0
    tmp_j = 0
    tmp_k = 0
    for n1 in range(1, a_r2+1):
        for n2 in range(n1+1, a_r2+1):
            if (      n1**2 +       n2**2 == a_r2):
                tmp_i += 1
            if (      n1**2 + 2.0 * n2**2 == a_r2 or
                2.0 * n1**2 +       n2**2 == a_r2):
                tmp_j += 1
            
            for n3 in range(n2+1, a_r2+1):
                if (n1**2 + n2**2 + n3**2 == a_r2):
                    tmp_k += 1
    
    ret += tmp_i * 12 * 2
    ret += tmp_j *  8 * 3
    ret += tmp_k *  8 * 6
    
    return ret

def num_point_r2_rev2(a_r2):
    """
    Return the number of r^2 which satisfies x^2 + y^2 + z^2 = r^2, where x, y, z are in Z.
    ~~~ Modified version2 ~~~
    """
    if (a_r2 == 0):
        return 1
    
    ret = 0
    
    if(int(sqrt(a_r2/1.0))**2 == a_r2/1.0):
        ret += 6
    if(int(sqrt(a_r2/2.0))**2 == a_r2/2.0):
        ret += 12
    if(int(sqrt(a_r2/3.0))**2 == a_r2/3.0):
        ret += 8
    
    tmp_i = 0
    tmp_j = 0
    tmp_k = 0
    for n1 in range(2, int(ceil(sqrt(a_r2)))):
        n2 = sqrt(a_r2 - n1**2)
        if (int(n2) == n2 and n1 > n2):
            tmp_i += 1
        
        n2 = sqrt((a_r2 - n1**2)/2.0)
        if (int(n2) == n2 and n1 > n2):
            tmp_j += 1
        if (a_r2 > 2.0*n1**2):
            n2 = sqrt(a_r2 - 2.0*n1**2)
            if (int(n2) == n2 and n1 > n2):
                tmp_j += 1
        
        for n2 in range(3, int(ceil(sqrt(a_r2 - n1**2)))):
            n3 = sqrt(a_r2-n1**2-n2**2)
            if (int(n3) == n3 and n2 > n1 > n3):
                tmp_k += 1
    
    return ret + (tmp_i * 12 * 2) + (tmp_j * 8 * 3) + (tmp_k * 8 * 6)

def zeta_1_yamazaki(a_q_sqr, precision = 10e-9):
    """
    The function to calculate the 'sum_{n in Z^3} (n^2 - a_q_sqr)^{-1}'
    The method of calculation is refered from T. Yamazaki et al., Phys. Rev. D 70 (2004) 074513
    
    return: 'sum_{n in Z^3} (n^2 - a_q_sqr)^{-1}'
    """
    term1 = sum_1_Z3 (lambda a_n2: exp(-(a_n2-a_q_sqr)) / (a_n2-a_q_sqr), precision)
    term2 = quad     (lambda a_t: exp(a_t*a_q_sqr) * (pi/a_t)**1.5 * 
                      sum_1_Z3(lambda a_n2: exp(-pi**2 * a_n2 / a_t), precision), 0, 1)[0]
    term3 = sum_0_inf(lambda a_l: (pi**1.5 * a_q_sqr**a_l) / ((a_l-0.5)*factorial(a_l)), precision)
    
    return term1 + term2 + term3

def zeta_1(a_q_sqr, precision = 10e-9):
    """
    The function to calculate the 'sum_{n in Z^3} (n^2 - a_q_sqr)^{-1}'
    
    return: 'sum_{n in Z^3} (n^2 - a_q_sqr)^{-1}'
    """
    return sum_1_Z3(lambda a_n2: 1.0/(a_n2 - a_q_sqr), precision) - 1.0/a_q_sqr

def k_cot_d_Luscher(a_k_sqr, a_L, precision = 10e-9):
    """
    The function to calculate 'k cot[delta(k)]' by Luscher's finite volume formula
    
    return: k cot[delta(k)]
    
    Note: a_k_sqr and a_L should be in Lattice Unit.
    """
    if (a_k_sqr == 0.0):
        return 0.0
    
    l_q_sqr = (a_L / (2.0 * pi))**2 * a_k_sqr
    
    #return zeta_1(l_q_sqr, precision) / (a_L * pi)
    return zeta_1_yamazaki(l_q_sqr, precision) / (a_L * pi)
