# -*- coding: utf-8 -*-

"""
The module for the Misner's method.
Ref: 'Spherical harmonic decomposition on a cubic grid', 
...  Charles W Misner, Classical and Quantum Gravity 21 (2004) S243-S247.
"""

from numpy                      import sqrt, where, array, empty, conj, dot, sum
from numpy.linalg               import inv
from scipy.special              import eval_legendre 
from misc_QM.special_functions  import sph_harm_xyz
from lattice.make_r_coordinates import xyzdata_to_rdata

def cmpl_set_Rn(a_n, a_R, a_del, a_x, a_y, a_z):
    """
    The function for complete set in the radial direction.
    ~~~ Legendre polynomial version ~~~
    
    return: Rn(x,y,z; R,Delta)
    """
    l_r = sqrt(a_x**2 + a_y**2 + a_z**2)
    return (sqrt((2*a_n+1)/(2*a_del)) * 
            eval_legendre(a_n, (l_r-a_R)/a_del) / l_r)

def cmpl_set_Rn_zero(a_n, a_R, a_del):
    """
    The function for complete set in the radial direction (at r = R).
    ~~~ Legendre polynomial version ~~~
    
    return: Rn(R; R,Delta)
    """
    return (sqrt((2*a_n+1)/(2*a_del)) * 
            eval_legendre(a_n, 0.0) / a_R)

def cmpl_set_Ynlm(a_n, a_l, a_m, a_R, a_del, a_x, a_y, a_z):
    """
    The function for complete set in the 2-sphere.
    ~~~ Legendre polynomial version ~~~
    
    return: Ynlm(x,y,z; R,Delta)
    """
    return (cmpl_set_Rn (a_n, a_R, a_del, a_x, a_y, a_z) *
            sph_harm_xyz(a_l, a_m,        a_x, a_y, a_z))

def mesure_cubed_sphere(a_R, a_del, a_x, a_y, a_z):
    """
    The function to calculate the measure in the 2-cubed sphere.
    """
    l_r  = sqrt(a_x**2 + a_y**2 + a_z**2)
    l_dr = abs(l_r-a_R)
    if   (l_dr > a_del + 0.5):
        return 0.0
    elif (l_dr < a_del - 0.5):
        return 1.0
    else:
        return a_del + 0.5 - l_dr

class Sph_harm_deco_misner:
    def __init__(self, a_maxN, a_maxL, a_delta):
        self.maxN  = a_maxN
        self.maxL  = a_maxL
        self.delta = a_delta
        self.nlm_t = [(n,l,m) 
                      for n in range(self.maxN+1) 
                      for l in range(self.maxL+1) 
                      for m in range(-l,l+1)]
        self.N_nlm = len(self.nlm_t)
        self.i_nlm = {}
        for i in range(self.N_nlm):
            self.i_nlm[self.nlm_t[i]] = i
        
        self.have_data = False
    
    def set_data3D(self, a_data3D):
        self.data      = xyzdata_to_rdata(a_data3D, False)
        self.have_data = True
    
    def set_data_r(self, a_data_r):
        self.data      = a_data_r
        self.have_data = True
    
    def calc_coef(self, a_R):
        if (not self.have_data):
            print("\nERROR: The data have not been setted.\n")
            return None
        
        idcs = where((a_R - self.delta - sqrt(3)/2 <= self.data[0].real) & 
                     (a_R + self.delta + sqrt(3)/2 >  self.data[0].real))[0]
        Np   = len(idcs)
        if (Np == 0):
            return None
        
        P_r = self.data[0,idcs].real
        Psi = self.data[1,idcs]
        P_z = self.data[2,idcs].real
        P_y = self.data[3,idcs].real
        P_x = self.data[4,idcs].real
        
        Ynlm = array([cmpl_set_Ynlm(*self.nlm_t[inlm], a_R,self.delta, P_x,P_y,P_z) 
                      for inlm in range(self.N_nlm)])
        
        omg  = array([mesure_cubed_sphere(a_R,self.delta, P_x[i],P_y[i],P_z[i]) 
                      for i in range(Np)])
        
        Gmat = empty((self.N_nlm, self.N_nlm), dtype=complex)
        for inlm in range(self.N_nlm):
            for jnlm in range(inlm, self.N_nlm):
                Gmat[inlm,jnlm] = sum(conj(Ynlm[inlm]) * Ynlm[jnlm] * omg)
                if (inlm != jnlm):
                    Gmat[jnlm,inlm] = conj(Gmat[inlm,jnlm])
        
        Gmat_inv = inv(Gmat)
        adj_Ynlm = array([dot(Ynlm[:,i], Gmat_inv) for i in range(Np)])
        
        return array([sum(conj(adj_Ynlm[:,inlm])*Psi*omg) for inlm in range(self.N_nlm)])
    
    def get_sph_amp(self, a_coef, a_R, a_l, a_m):
        if (a_l > self.maxL):
            print("\nERROR: a_l = %d is too large than maxL = %d.\n" % (a_l, self.maxL))
            return None
        if (a_m < -a_l or a_l < a_m):
            print("\nERROR: Invalid value: a_m = %d.\n" % (a_m))
            return None
        
        return sum(array([a_coef[self.i_nlm[n,a_l,a_m]] * 
                          cmpl_set_Rn_zero(n, a_R,self.delta) for n in range(self.maxN+1)]))
    
    def calc_sph_amp_lm(self, a_coef, a_R, a_l, a_m):
        return self.get_sph_amp(self.calc_coef(a_R), a_R, a_l, a_m)
