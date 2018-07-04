# -*- coding: utf-8 -*-

"""The module to calculate phase shift from S-matrix."""

def calc_phase_Sii(a_Sii, error_crit = 1e-10, shift_disc = 45):
    """
    The function to calculate phase shift from S-matrix.
    ~~~ For diagonal elements of S-matrix ~~~
    
    For arguments,
    - a_Sii[#.energy] (1-dim ndarray, dtype=complex)
    
    return: np.array([r_phase[#.energy], r_inelasticity[#.energy]])
    
    Note1: #.energy is got from a_Sii.
    Note2: The phase shift is defined as 'log(Sii/|Sii|).imag/2.0'
    .....: in the unit of degrees, where the inelasticity is |Sii|.
    Note3: If log(Sii/|Sii|).real/2.0 > error_criterion, return (None, None)
    """
    from cmath import log
    from numpy import array, empty, pi, sign
    
    r_inela = array([abs(a_Sii[iE]            ) for iE in range(len(a_Sii))])
    l_phase = array([log(a_Sii[iE]/r_inela[iE]) for iE in range(len(a_Sii))])
    r_phase = empty(len(a_Sii))
    
    if (l_phase[0].real*90.0/pi > error_crit):
        print("\nERROR: log(Sii/|Sii|).real/2.0 > %e\n" % criterion)
        return (None, None)
    
    phs_tmp0   = l_phase[0].imag * 90.0 / pi
    r_phase[0] = phs_tmp0
    shift      = 0.0
    
    for iE in range(1, len(a_Sii)):
        if (l_phase[iE].real*90.0/pi > error_crit):
            print("\nERROR: log(Sii/|Sii|).real/2.0 > %e\n" % criterion)
            return (None, None)
        
        phs_tmp1 = l_phase[iE].imag * 90.0 / pi
        
        if (sign(phs_tmp0) != sign(phs_tmp1)):
            if (abs(phs_tmp0-phs_tmp1) > shift_disc):
                shift += -sign(phs_tmp1) * 180
        
        r_phase[iE] = phs_tmp1 + shift
        phs_tmp0    = phs_tmp1
    
    return array([r_phase, r_inela])

def within_one(a_num):
    """
    Sub-function to avoid the error about the asin or acos.
    
    0 <= abs(a_num) <= 1 --> return     a_num  ( 0 ~~ 1)
    1 <  abs(a_num) <  2 --> return int(a_num) (-1 or 1)
    2 <= abs(a_num)      --> print the error messages
    """
    
    if   (     abs(a_num) <= 1):
        return     a_num
    elif (1 <  abs(a_num) <  2):
        return int(a_num)
    else:
        print("\nERROR: In acos(x) or asin(x), x >= 2 happened.\n")
        return None

############################################
### The functions bellow will be deleted ###

from math import sqrt, sin, cos, asin, acos, atan, pi

def calc_phase_Nch1(a_Smat):
    """
    The function to calculate phase shift from S-matrix.
    ~~~ For #.channel == 1 ~~~
    
    For arguments,
    - a_Smat (scalar)
    
    return: phase_shift
    
    Note1:      a_Smat is supposed to complex-type.
    Note2: phase_shift is supposed to   float-type.
    """
    
    phs_s = asin(a_Smat.imag) * 90 / pi   # <= Works well for unbound
    phs_c = acos(a_Smat.real) * 90 / pi   # <= Works well for   bound
    
    return (phs_s, abs(a_Smat))

def calc_phase_Nch2(a_Smat):
    """
    The function to calculate phase shift from S-matrix.
    ~~~ For #.channel == 2 ~~~
    
    For arguments,
    - a_Smat[2, 2] (2-dim ndarray)
    
    return: (phase_shift[2], mixing-angle, inelasticity)
    
    Note1:      a_Smat is supposed to complex-type.
    Note2: phase_shift is supposed to   float-type.
    """
    
    #Tan2theta   = -(a_Smat[0,1]*a_Smat[1,0]) / (a_Smat[0,0]*a_Smat[1,1]) # Not necessary
    Cos2theta00 = abs(a_Smat[0,0])
    Cos2theta11 = abs(a_Smat[1,1])
    Sin2theta01 = abs(a_Smat[0,1])
    Sin2theta10 = abs(a_Smat[1,0])
    
    phs1_tmp00 = a_Smat[0,0] / Cos2theta00
    
    phs1_s = asin(within_one(phs1_tmp00.imag)) / 2.0   # <= Works well for unbound
    phs1_c = acos(within_one(phs1_tmp00.real)) / 2.0   # <= Works well for   bound
        
    if (abs(1-Cos2theta00) < 1e-6): # below the threshold
        r_phase2 = 0.0
        r_mixAng = 0.0
    else:
        phs2_tmp11 = a_Smat[1,1] / Cos2theta11
        
        phs2_s = asin(within_one(phs2_tmp11.imag)) / 2.0   # <= Works well for unbound
        phs2_c = acos(within_one(phs2_tmp11.real)) / 2.0   # <= Works well for   bound
        
        #phs12_tmp1 = a_Smat[0,1] / Sin2theta01   # <=             For Debug
        #phs12_tmp2 = a_Smat[1,0] / Sin2theta10   # <=             For Debug
        #phs12_s_1  = asin(-phs12_tmp1.real)      # <= Works well, For Debug
        #phs12_s_2  = asin(-phs12_tmp2.real)      # <= Works well, For Debug 
        #phs12_c_1  = acos( phs12_tmp1.imag)      # <= Works well, For Debug
        #phs12_c_2  = acos( phs12_tmp2.imag)      # <= Works well, For Debug
        
        mixA_s_1 = asin(within_one(Sin2theta01)) / 2.0   # <= Works Well, but not convenience
        mixA_s_2 = asin(within_one(Sin2theta10)) / 2.0   # <= Works Well, but not convenience
        mixA_c_1 = acos(within_one(Cos2theta00)) / 2.0   # <= Works well, but not convenience
        mixA_c_2 = acos(within_one(Cos2theta11)) / 2.0   # <= Works well, but not convenience
        
        if (abs(sin(phs1_s+phs2_s)) > 1.0/sqrt(2.0)):
            mixA = asin(within_one(a_Smat[0,1].real / -sin(phs1_s+phs2_s))) / 2.0   # <= Works very well
        else:
            mixA = asin(within_one(a_Smat[0,1].imag /  cos(phs1_s+phs2_s))) / 2.0   # <= Works very well
        
        r_phase2 = phs2_s * 180 / pi
        r_mixAng = mixA   * 180 / pi
    
    r_phase1 = phs1_s * 180 / pi
    r_inela  = Cos2theta00
    
    return ((r_phase1, r_phase2), r_mixAng, r_inela)
