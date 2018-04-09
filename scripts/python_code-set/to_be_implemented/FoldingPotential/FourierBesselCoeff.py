# -*- coding: utf-8 -*-

"""The module for Definition of the Fourier-Bessel Coefficients"""

from numpy import array, pi, sin, cos, sum, sqrt

def params_FBcoe(a_NA):
    """
    The function for the definition of parameters of the Fourier-Bessel Coefficients.
    
    return: array([R_cut, a_1, a_2, ...])
    
    Note: The definition of parameters is referred from
    .     'Nuclear charge-density-distribution parameters 
    .      from elastic electron scattering'
    .     H. De Vries, C. W. De Jager and C. De Vries
    .     Atomic Data and Nuclear Data Tables, 36, 3 (1987) p.495-536
    """
    
    if   (a_NA ==  12): # For  12 C , RMS = 2.472(15) fm
        RetArray = array([8.0,
                          0.15721e-1, 0.38732e-1, 0.36808e-1, 0.14671e-1, -0.43277e-2,
                          -0.97752e-2, -0.68908e-2, -0.27631e-2, -0.63568e-3, 0.71809e-4,
                          0.18441e-3, 0.75066e-4, 0.51069e-4, 0.14308e-4, 0.23170e-5,
                          0.68465e-6])
        
    elif (a_NA ==  28): # For  28 Si, RMS = 3.085(17) fm
        RetArray = array([8.0,
                          0.33495e-1, 0.59533e-1, 0.20979e-1, -0.16900e-1, -0.14998e-1,
                          -0.93248e-3, 0.33266e-2, 0.59244e-3, -0.40013e-3, 0.12242e-3,
                          -0.12994e-4, -0.92784e-5, 0.72595e-5, -0.42096e-5])
        
    elif (a_NA ==  40): # For  40 Ca, RMS = 3.450(10) fm
        RetArray = array([8.0,
                          0.44846e-1, 0.61326e-1, -0.16818e-2, -0.26217e-1, -0.29725e-2,
                          0.85534e-2, 0.35322e-2, -0.48258e-3, -0.39346e-3, 0.20338e-3,
                          0.25461e-4, -0.17794e-4, 0.67394e-5, -0.21033e-5])
        
    elif (a_NA ==  58): # For  58 Ni, RMS = 3.769(13) fm
        RetArray = array([9.0,
                          0.44880e-1, 0.64756e-1, -0.27899e-2, -0.37016e-1, -0.71915e-2,
                          0.13594e-1, 0.66331e-2, -0.14095e-2, -0.10141e-2, 0.38616e-3,
                          -0.13871e-3, 0.47788e-4, -0.15295e-4, 0.59131e-5, -0.67880e-5])

    elif (a_NA ==  90): # For  90 Zr, RMS = 4.258( 8) fm
        RetArray = array([10.0,
                          0.46188e-1, 0.61795e-1, -0.12315e-1, -0.36915e-1, 0.25175e-2,
                          0.15234e-1, -0.55146e-3, -0.60631e-2, -0.12198e-2, 0.36200e-3,
                          -0.16466e-3, 0.53305e-4, -0.50873e-5, -0.85658e-5, 0.86095e-5])
        
    elif (a_NA == 208): # For 208 Pb, RMS = 5.499( 1) fm
        RetArray = array([11.0,
                          0.62732e-1, 0.38542e-1, -0.55105e-1, -0.26990e-2, 0.31016e-1,
                          -0.99486e-2, -0.93012e-2, 0.76653e-2, 0.20885e-2, -0.17840e-2,
                          0.74876e-4, 0.32278e-3, -0.11353e-3])
        
    else:
        print("\nWARNING: #.nucleus = %d has not been implemented yet.\n" % a_NA)
        RetArray = None
    
    return RetArray

def norm_FBcoe(a_NA):
    """
    The function for the Normalization of the Fourier-Bessel Coefficients.
    
    return: Int d^3 r rho_FB(r)
    """
    
    Rcut = params_FBcoe(a_NA)[0]
    a_nu = params_FBcoe(a_NA)[1:]
    nu   = array([i for i in range(1, len(a_nu)+1)])
    
    return 4.0 * Rcut**3 / pi**2 * sum(a_nu/nu**3 * (sin(nu*pi) - nu*pi * cos(nu*pi)))

def rms_FBcoe(a_NA):
    """
    The function for the Root-Mean-Square of the Fourier-Bessel Coefficients.
    
    return: Sqrt( < r^2 > )
    """
    
    Rcut = params_FBcoe(a_NA)[0]
    a_nu = params_FBcoe(a_NA)[1:]
    nu   = array([i for i in range(1, len(a_nu)+1)])
    
    Rsq  = (4.0 * Rcut**5 / pi**4 * 
            sum(a_nu/nu**5 * (3.0   * ((nu*pi)**2 - 2.0) * sin(nu*pi) - 
                              nu*pi * ((nu*pi)**2 - 6.0) * cos(nu*pi))))
    
    return sqrt(Rsq/norm_FBcoe(a_NA))
