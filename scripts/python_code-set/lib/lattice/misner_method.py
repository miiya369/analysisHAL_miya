# -*- coding: utf-8 -*-

"""
The module for the Misner's method.
Ref: 'Spherical harmonic decomposition on a cubic grid', 
...  Charles W Misner, Classical and Quantum Gravity 21 (2004) S243-S247.
"""

from numpy                      import sqrt, where, array, empty, conj, dot, sum
from numpy.linalg               import inv
from scipy.special              import eval_legendre, legendre
from multiprocessing            import Pool

from misc_QM.special_functions  import sph_harm_xyz

class Sph_harm_deco_misner:
    def __init__(self, a_maxN, a_maxL, a_delta):
        """
        Brief: Constructor of the class.
        
        For arguments,
        - a_maxN  (scalar integer): #.basis for the complete set in the radial direction
        .                           to be used for the calculation.
        
        - a_maxL  (scalar integer): #.basis for the spherical harmonics to be used,
        .                           where M is selected as -L <= M <= L for each L.
        
        - a_delta (scalar  float ): Thickness of the spherical shell (S), defined as
        .                           S = { (r,theta,phi) | R - delta <= r < R + delta } 
        .                           for each R.
        """
        self.maxN     = a_maxN.real
        self.maxL     = a_maxL.real
        self.delta    = a_delta.real
        self.calcR    = None
        self.oRlm_all = None
        
        # Index table for (n,l,m)
        self.nlm_t = [(n,l,m) 
                      for n in range(self.maxN+1) 
                      for l in range(self.maxL+1) 
                      for m in range(-l,l+1)]
        self.N_nlm = len(self.nlm_t)
        self.i_nlm = {}
        for i in range(self.N_nlm):
            self.i_nlm[self.nlm_t[i]] = i
        
        # Index table for (l,m)
        self.lm_t = [(l,m) 
                      for l in range(self.maxL+1) 
                      for m in range(-l,l+1)]
        self.N_lm = len(self.lm_t)
        self.i_lm = {}
        for i in range(self.N_lm):
            self.i_lm[self.lm_t[i]] = i
        
        self.cmpl_set_Rn          = self.cmpl_set_Rn_Legendre
        self.cmpl_set_Rn_at_R     = self.cmpl_set_Rn_Legendre_at_R
        self.cmpl_set_lap_Rn_at_R = self.cmpl_set_lap_Rn_Legendre_at_R
        self.measure_cubed_sphere = self.measure_cubed_sphere_simple
    
    def set_Rn(self, a_type):
        """
        Brief: Set the functions for the complete set in r-direction.
        """
        if (a_type == 0):
            print("Set the complete set in r-direction: Legendre polynomial type")
            self.cmpl_set_Rn          = self.cmpl_set_Rn_Legendre
            self.cmpl_set_Rn_at_R     = self.cmpl_set_Rn_Legendre_at_R
            self.cmpl_set_lap_Rn_at_R = self.cmpl_set_lap_Rn_Legendre_at_R
        else:
            print("The type = %d has not been implemented." % a_type)
    
    def set_measure(self, a_type):
        """
        Brief: Set the function for the measure in the cubed sphere.
        """
        if (a_type == 0):
            print("Set the measure in the cubed sphere: simple type")
            self.measure_cubed_sphere = self.measure_cubed_sphere_simple
        else:
            print("The type = %d has not been implemented." % a_type)
    
    def set_coords(self, a_R, a_rzyx_coords, a_calc_orthogonal_basis = True):
        """
        Brief: Set data points (coordinats) to be calculated.
        
        For arguments,
        - a_R          [   #.output-data] (1-dim ndarray):
        . Values of the r-coordinate for output data.
        
        - a_rzyx_coords[4, #.data-points] (2-dim ndarray):
        . Values of the r,z,y,x-coordinates for input data.
        . This argument must have the following structure,
        . 
        . a_rzyx_coords[0, i] -> r-coordinate for i-th data point,
        . a_rzyx_coords[1, i] -> z-coordinate for i-th data point,
        . a_rzyx_coords[2, i] -> y-coordinate for i-th data point,
        . a_rzyx_coords[3, i] -> x-coordinate for i-th data point,
        . 
        . for i = 0 ~ N-1, where N is the number of data points.
        
        - a_calc_orthogonal_basis (boolean): True for calculation of orthogonal basis
        """
        self.calcR  = a_R.real
        self.Np_out = len(self.calcR)
        self.irzyx  = a_rzyx_coords.real
        
        r_margin   = sqrt(3)/2
        self.idcs  = [where((self.calcR[i] - self.delta - r_margin <= self.irzyx[0]) & 
                            (self.calcR[i] + self.delta + r_margin >  self.irzyx[0]))[0] 
                      for i in range(self.Np_out)]
        self.Np_in = array([len(self.idcs[i]) for i in range(self.Np_out)])
        
        if (a_calc_orthogonal_basis):
            self.calc_orthogonal_basis()
    
    def calc_orthogonal_basis_in_Si(self, i):
        """
        Brief: Calculate orthogonal basis in the spherical shell 'S_i' defined as
        .      S_i = { (r,theta,phi) | R_i - delta <= r < R_i + delta },
        .      where R_i is i-th r-coordinate for the output data.
        
        return: oRlm[ 2, self.N_lm, self.Np_in[i] ] (3-dim ndarray, dtype=complex)
        .       oRlm[ 0 ] --> omega * Sum_n [              Rn (r=R) * conj(Y^{nlm}) ]
        .       oRlm[ 1 ] --> omega * Sum_n [ Laplacian of Rn (r=R) * conj(Y^{nlm}) ],
        .       where Y^{nlm} is the adjoint basis defined in Misner's paper.
        """
        if (self.calcR is None):
            print("\nERROR: The data points (coordinats) have not been setted.\n")
            return None
        if (i >= self.Np_out):
            print("\nERROR: Index overflow.\n")
            return None
        
        crd_r = self.irzyx[0,self.idcs[i]]
        crd_z = self.irzyx[1,self.idcs[i]]
        crd_y = self.irzyx[2,self.idcs[i]]
        crd_x = self.irzyx[3,self.idcs[i]]
        
        Ynlm = array([self.cmpl_set_Ynlm(*self.nlm_t[inlm], self.calcR[i], 
                                          crd_r,crd_x,crd_y,crd_z) 
                      for inlm in range(self.N_nlm)])
        
        omg  = array([self.measure_cubed_sphere(self.calcR[i], 
                                                crd_x[ip],crd_y[ip],crd_z[ip])
                      for ip in range(self.Np_in[i])])
        
        Gmat = empty((self.N_nlm, self.N_nlm), dtype=complex)
        for inlm in range(self.N_nlm):
            for jnlm in range(inlm, self.N_nlm):
                Gmat[inlm,jnlm] = sum(conj(Ynlm[inlm]) * Ynlm[jnlm] * omg)
                if (inlm != jnlm):
                    Gmat[jnlm,inlm] = conj(Gmat[inlm,jnlm])
        
        conj_adj_Ynlm = conj(dot(Ynlm.T, inv(Gmat)))
        
        oRlm = array([sum(array([self.cmpl_set_Rn_at_R(n, self.calcR[i]) *
                                 conj_adj_Ynlm[:,self.i_nlm[n,l,m]]
                                 for n in range(self.maxN)]), axis=0) * omg
                      for l in range(self.maxL+1) for m in range(-l,l+1)])
        
        lap_oRlm = array([sum(array([self.cmpl_set_lap_Rn_at_R(n, self.calcR[i]) *
                                     conj_adj_Ynlm[:,self.i_nlm[n,l,m]]
                                     for n in range(self.maxN)]), axis=0) * omg
                          for l in range(self.maxL+1) for m in range(-l,l+1)])
        
        return array([oRlm, lap_oRlm])
    
    def calc_orthogonal_basis(self, a_Nprocs = 1):
        """
        Brief: Calculate orthogonal basis in all spherical shells for the output data.
        
        For arguments,
        - a_Nprocs (scalar integer): #.processes used for the calculation.
        """
        if (self.calcR is None):
            print("\nERROR: The data points (coordinats) have not been setted.\n")
            return None
        
        if (a_Nprocs == 1):
            self.oRlm_all = [self.calc_orthogonal_basis_in_Si(i) 
                             for i in range(self.Np_out)]
        else:
            with Pool(a_Nprocs) as proc:
                self.oRlm_all = proc.map(self.calc_orthogonal_basis_in_Si, 
                                         [i for i in range(self.Np_out)])
    
    def get_sph_amp(self, a_rdata, a_l, a_m, a_get_laplacian = False):
        """
        Brief: Get an amplitude of the spherical harmonics for given data.
        
        For arguments,
        - a_rdata[#.data-points] ( 1-dim ndarray): Input data
        - a_l, a_m               (scalar integer): (l,m)-set to be returned.
        - a_get_laplacian        (       boolean): True for return the laplacian.
        
        Note: The order of the elements of a_rdata must be the same with that of
        .     a_rzyx_coords which are arguments of self.set_coords.
        .     That is, when the value of i-th data is a_rdata[i], the r-coordinate of
        .     this data shold be a_rzyx_coords[0, i], the z-coordinate of this data 
        .     shold be a_rzyx_coords[0, i], and so on.
        
        return: amplitude for l,m [#.output-data] (1-dim ndarray, dtype=complex).
        """
        if (self.oRlm_all is None):
            print("\nERROR: The orthogonal basis have not been calculated.\n")
            return None
        if (a_l > self.maxL):
            print("\nERROR: a_l = %d is too large than maxL = %d.\n" % (a_l, self.maxL))
            return None
        if (a_m < -a_l or a_l < a_m):
            print("\nERROR: Invalid value: a_m = %d.\n" % a_m)
            return None
        
        if (a_get_laplacian):
            return array([sum(self.oRlm_all[i][1, self.i_lm[a_l,a_m]] * 
                              a_rdata[self.idcs[i]]) for i in range(self.Np_out)])
        else:
            return array([sum(self.oRlm_all[i][0, self.i_lm[a_l,a_m]] * 
                              a_rdata[self.idcs[i]]) for i in range(self.Np_out)])
    
    def cmpl_set_Rn_Legendre(self, a_n, a_R, a_r):
        """
        Brief: The function for complete set in the radial direction.
        .      ~~~ Legendre polynomial version ~~~
        
        return: Rn(r; R,Delta).
        """
        return (sqrt((2*a_n+1)/(2*self.delta)) * 
                eval_legendre(a_n, (a_r-a_R)/self.delta) / a_r)
    
    def cmpl_set_Rn_Legendre_at_R(self, a_n, a_R):
        """
        Brief: The function for complete set in the radial direction (at r = R).
        .      ~~~ Legendre polynomial version ~~~
        
        return: Rn(R; R,Delta).
        """
        return sqrt((2*a_n+1)/(2*self.delta)) * eval_legendre(a_n, 0.0) / a_R
    
    def cmpl_set_lap_Rn_Legendre_at_R(self, a_n, a_R):
        """
        Brief: laplacian derivative of orthogonal basis function in radial direction (at r = R).
        .      ~~~ Legendre polynomial version ~~~
        
        return: 1/r d^2/dr^2 r Rn(r) at r=R.
        """
        return (sqrt((2*a_n+1)/(2*self.delta)) * 
                legendre(a_n).deriv(m=2)(0.0) / (a_R * self.delta**2))
    
    def cmpl_set_Ynlm(self, a_n, a_l, a_m, a_R, a_r, a_x, a_y, a_z):
        """
        Brief: The function for complete set in the 2-sphere.
        
        return: Ynlm(x,y,z; R,Delta).
        """
        return (self.cmpl_set_Rn(a_n, a_R, a_r) *
                sph_harm_xyz    (a_l, a_m, a_x, a_y, a_z))
    
    def measure_cubed_sphere_simple(self, a_R, a_x, a_y, a_z):
        """
        Brief: The function to calculate the measure in the 2-cubed spheres.
        .      ~~~ Simple version ~~~
        
        return: omega(x,y,z; R,Delta).
        """
        l_dr = abs(sqrt(a_x**2 + a_y**2 + a_z**2)-a_R)
        if   (l_dr > self.delta + 0.5):
            return 0.0
        elif (l_dr < self.delta - 0.5):
            return 1.0
        else:
            return self.delta + 0.5 - l_dr
