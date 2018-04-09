# -*- coding: utf-8 -*-

"""The module for Definition of the density functions. (w/o rho_0)"""

from numpy import vectorize, sum, array, exp, pi, sin

sph_j0  = lambda x: sin(x) / x if x != 0 else 1.0
Vsph_j0 = vectorize(sph_j0)

def dens_FBcoe(a_r, a_Params_FBcoe):
    """
    The function for the Fourier-Bessel Coefficients
    
    For arguments,
    - a_Params_FBcoe[R_cut, a_1, a_2, ...]
    
    return: Sum_{n=1}^N a_n * j_0(n pi r / R_cut)
    """
    
    Rcut = a_Params_FBcoe[0]
    a_nu = a_Params_FBcoe[1:]
    nu   = array([i for i in range(1, len(a_nu)+1)])
    
    if (a_r <= Rcut):
        return sum(a_nu * Vsph_j0(nu*pi * a_r / Rcut))
    else:
        return 0.0

def dens_gaussian(a_r, a_R_A):
    """
    The function for the Gaussian-type density (w/o rho_0).
    
    return: exp[-(a_r/R_A)^2]
    """
    
    return exp(-pow(a_r / a_R_A, 2))

def dens_woods_saxon(a_r, a_R_A, a_dif_A):
    """
    The function for the Woods-Saxon-type density (w/o rho_0).
    
    return: (1 + exp[(a_r-R_A)/dif_A])^(-1)
    
    Note: In practice, exp[(R_A-a_r)/dif_A] / (exp[(R_A-a_r)/dif_A] + 1)
          is returned for numerical stability.
    """
    
    return (exp ((a_R_A - a_r) / a_dif_A) / 
            (exp((a_R_A - a_r) / a_dif_A) + 1.0))

def calc_rho0(a_A, a_func_dens, a_Params):
    """
    The function to calculate the rho_0 for given density function.
    
    For arguments,
    - a_func_dens       (object of function)
    - a_Params[#.param] (       1-dim array)
    
    return: rho_0
    """
    
    l_flg = True # For Debug
    
    if (a_func_dens.__name__ == 'dens_woods_saxon' and l_flg):
        from mpmath import polylog
        r_rho0 = a_A / (-8*pi * a_Params[1]**3 * polylog(3, -exp(a_Params[0]/a_Params[1])))
        
        #r_rho0 = 3.0*a_A / (4*pi * (a_Params[0]**3 + (pi*a_Params[1])**2 * a_Params[0])) # Approximate velue
        
    elif (a_func_dens.__name__ == 'dens_gaussian' and l_flg):
        r_rho0 = a_A / (pow(pi, 1.5) * a_Params[0]**3)
        
    else:
        from scipy.integrate import quad
        from numpy           import inf
        Ret_rho0 = quad(lambda r: r**2 * a_func_dens(r, *a_Params), 0.0, inf)[0]
        r_rho0   = a_A / (4.0 * pi * Ret_rho0)
    
    return r_rho0
