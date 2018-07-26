# -*- coding: utf-8 -*-

"""The module to calculate the potential by HAL QCD method on the lattice."""

from lattice.calc_laplacian import lap, lap_4th_proc

def calc_potential_t2(a_Rm, a_R0, a_Rp, a_mass):
    """
    The function to calculate the potential by HAL QCD method on the lattice.
    ~~~ With Laplacian-term + 1st t-deriv. term + 2nd t-deriv. term ~~~
    
    For arguments,
    - a_Rm, a_R0, a_Rp[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - a_mass[2]                                   (1-dim ndarray)
    ( a_mass = [mass1, mass2] )
    
    return: r_pot[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    """
    l_mu = (a_mass[0] * a_mass[1]) / (a_mass[0] + a_mass[1])
    l_dl = (a_mass[0] - a_mass[1]) / (a_mass[0] + a_mass[1])
    
    return (lap_4th_proc(a_R0) / (2.0*l_mu) +
            (a_Rm - a_Rp) / 2.0 +
            (a_Rm + a_Rp - 2.0*a_R0) * (1.0+3.0*l_dl**2) / (8.0*l_mu)
            ) / a_R0

def calc_potential_t1(a_Rm, a_R0, a_Rp, a_mass):
    """
    The function to calculate the potential by HAL QCD method on the lattice.
    ~~~ With Laplacian-term + 1st t-deriv. term ~~~
    
    For arguments,
    - a_Rm, a_R0, a_Rp[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - a_mass[2]                                   (1-dim ndarray)
    ( a_mass = [mass1, mass2] )
    
    return: r_pot[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    """
    l_mu = (a_mass[0] * a_mass[1]) / (a_mass[0] + a_mass[1])
    l_dl = (a_mass[0] - a_mass[1]) / (a_mass[0] + a_mass[1])
    
    return (lap_4th_proc(a_R0) / (2.0*l_mu) +
            (a_Rm - a_Rp) / 2.0
            ) / a_R0

def calc_potential_t2_num(a_Rm, a_R0, a_Rp, a_mass):
    """
    The function to calculate the potential by HAL QCD method on the lattice.
    ~~~ With Laplacian-term + 1st t-deriv. term + 2nd t-deriv. term ~~~
    ~~~ Output only the numerator ~~~
    
    For arguments,
    - a_Rm, a_R0, a_Rp[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - a_mass[2]                                   (1-dim ndarray)
    ( a_mass = [mass1, mass2] )
    
    return: r_pot[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    """
    l_mu = (a_mass[0] * a_mass[1]) / (a_mass[0] + a_mass[1])
    l_dl = (a_mass[0] - a_mass[1]) / (a_mass[0] + a_mass[1])
    
    return (lap_4th_proc(a_R0) / (2.0*l_mu) +
            (a_Rm - a_Rp) / 2.0 +
            (a_Rm + a_Rp - 2.0*a_R0) * (1.0+3.0*l_dl**2) / (8.0*l_mu)
            )

def calc_potential_t1_num(a_Rm, a_R0, a_Rp, a_mass):
    """
    The function to calculate the potential by HAL QCD method on the lattice.
    ~~~ With Laplacian-term + 1st t-deriv. term ~~~
    ~~~ Output only the numerator ~~~
    
    For arguments,
    - a_Rm, a_R0, a_Rp[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - a_mass[2]                                   (1-dim ndarray)
    ( a_mass = [mass1, mass2] )
    
    return: r_pot[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    """
    l_mu = (a_mass[0] * a_mass[1]) / (a_mass[0] + a_mass[1])
    l_dl = (a_mass[0] - a_mass[1]) / (a_mass[0] + a_mass[1])
    
    return (lap_4th_proc(a_R0) / (2.0*l_mu) +
            (a_Rm - a_Rp) / 2.0
            )
