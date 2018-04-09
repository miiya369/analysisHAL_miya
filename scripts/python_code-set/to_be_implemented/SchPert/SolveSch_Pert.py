# -*- coding: utf-8 -*-

"""The module to solve the Schrodinger equation by perturbation."""

from scipy.integrate import quad
from numpy           import array, inf, exp, sqrt, pi

hbar_c = 197.327053

def solve_sch_0th_pert(a_Range, a_Coeff, mass, a_FuncV, a_ParamV):
    """
    The function to solve the Schrodinger equation by perturbation.
    ~~~ For Original (0-th) Energy ~~~
    
    For arguments,
    - a_Range [#.base]  (1-dim array)
    - a_Coeff [#.base]  (1-dim array)
    - a_ParamV[#.param] (1-dim array)
    - a_FuncV           (pointer to the potential function)
    
    return: (0-th Energy) as E_0 = <0|(H_0 + V)|0> / <0|0>
    
    Note: #.base is got from a_Range.
    """
    
    Nbase = len(a_Range)
    
    H_0 = 0.0; Nrm = 0.0; Pot = 0.0
    for i in range(Nbase):
        for j in range(i, Nbase):
            F = 1.0 if i == j else 2.0
            
            aij2 =         a_Range[i]**2 + a_Range[j]**2
            Cij  =         a_Coeff[i]    * a_Coeff[j]
            Aij2 = aij2 / (a_Range[i]**2 * a_Range[j]**2)
            
            H_0 += F * Cij / (aij2 * Aij2**1.5)
            Nrm += F * Cij /  Aij2**1.5
            Pot += F * Cij * quad(lambda r: r**2 * a_FuncV(r, *a_ParamV) * exp(-Aij2*r**2), 
                                  0, inf)[0]
    
    return (3.0*hbar_c**2/mass * H_0 + 4.0/sqrt(pi) * Pot) / Nrm

def solve_sch_1st_pert(a_Range, a_Coeff, a_FuncV, a_ParamV):
    """
    The function to solve the Schrodinger equation by perturbation.
    ~~~ For 1st-order Energy of the perturbation ~~~
    
    For arguments,
    - a_Range [#.base]  (1-dim array)
    - a_Coeff [#.base]  (1-dim array)
    - a_ParamV[#.param] (1-dim array)
    - a_FuncV           (pointer to the potential function)
    
    return: (1st Energy) as E_1 = <0|V'|0> / <0|0>
    
    Note: #.base is got from a_Range.
    """
    
    Nbase = len(a_Range)
    
    Nrm = 0.0; Pot = 0.0
    for i in range(Nbase):
        for j in range(i, Nbase):
            F = 1.0 if i == j else 2.0
            
            aij2 =         a_Range[i]**2 + a_Range[j]**2
            Cij  =         a_Coeff[i]    * a_Coeff[j]
            Aij2 = aij2 / (a_Range[i]**2 * a_Range[j]**2)
            
            Nrm += F * Cij / Aij2**1.5
            Pot += F * Cij * quad(lambda r: r**2 * a_FuncV(r, *a_ParamV) * exp(-Aij2*r**2), 
                                  0, inf)[0]
    
    return 4.0/sqrt(pi) * Pot / Nrm

def solve_sch_1st_pert_Coulomb_FBcoe(a_Range, a_Coeff, a_NA):
    """
    The function to solve the Schrodinger equation by perturbation.
    ~~~             For 1st-order Energy of the perturbation              ~~~
    ~~~ For the folding potential of Coulomb + Fourier-Bessel Coefficient ~~~
    
    For arguments,
    - a_Range [#.base]  (1-dim array)
    - a_Coeff [#.base]  (1-dim array)
    
    return: (1st Energy) as E_1 = <0|Vc|0> / <0|0>
    
    Note: #.base is got from a_Range.
    """
    
    from FoldingPotential.FourierBesselCoeff import params_FBcoe
    from FoldingPotential.PotentialFunc      import fpot_Coulomb_dens_FBcoe
    
    Nbase = len(a_Range)
    
    Nrm = 0.0; Pot = 0.0
    for i in range(Nbase):
        for j in range(i, Nbase):
            F = 1.0 if i == j else 2.0
            
            aij2 =         a_Range[i]**2 + a_Range[j]**2
            Cij  =         a_Coeff[i]    * a_Coeff[j]
            Aij2 = aij2 / (a_Range[i]**2 * a_Range[j]**2)
            
            Nrm += F * Cij / Aij2**1.5
            Pot += F * Cij * quad(lambda r: r**2 * exp(-Aij2*r**2) *
                                  fpot_Coulomb_dens_FBcoe(r, params_FBcoe(a_NA)),
                                  0, inf)[0]
    
    return 4.0/sqrt(pi) * Pot / Nrm
