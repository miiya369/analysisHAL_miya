# -*- coding: utf-8 -*-

"""
The module to solve the Schrodinger equation by using Gauss expansion method
with folding-potential (Woods-Saxon type).
"""

from numpy           import array
from numpy.linalg    import norm
from scipy.linalg    import eigh
from multiprocessing import Pool

from sch_gauss_exp.set_matrix_GEM_fpot_WS import *

def fproc_wrapper_make_Ham_mat(args):
    a_iconf       = args[0]
    a_range       = args[1]
    a_params      = args[2]
    a_mass        = args[3]
    a_args_dens   = args[4]
    a_func_name   = args[5]
    a_verbose_flg = args[6]
    
    if (re.match('^[1-9]G$', a_func_name) is not None):
        r_Ham = set_Ham_mat_gauss_fpot_WS(a_range, a_params, a_mass, a_args_dens)
    else:
        r_Ham = set_Ham_mat_fpot_WS(a_range, a_params, a_func_name, a_mass, a_args_dens)
    
    if (a_verbose_flg):
        print("# Make Hamiltonian matrixes for GEM... %3d end" % a_iconf)
    return r_Ham

def fproc_wrapper_solve_sch_GEM(args):
    a_iconf       = args[0]
    a_Ham_C       = args[1]
    a_Phi         = args[2]
    a_eval_only   = args[3]
    a_Nstat       = args[4]
    a_verbose_flg = args[5]
    
    r_eigens = eigh(a_Ham_C, a_Phi, eigvals_only=a_eval_only, eigvals=(0, a_Nstat-1))
    
    if (a_verbose_flg):
        print("# Solve Schrodinger equation by GEM... %3d end" % a_iconf)
    return r_eigens

def solve_sch_GEM_fpot_WS(a_Nstat, a_range, a_params, a_mass, a_args_dens, a_func_name, a_Np = 0, a_eval_only = True):
    """
    The function to solve the Schrodinger equation by using Gauss expansion method.
    
    For arguments,
    - a_range [#.base]  (1-dim ndarray)
    - a_params[#.param] (1-dim ndarray)
    - a_args_dens[4]    (1-dim ndarray)
    
    return: (eval, evec)
    - eval[#.state]         (1-dim ndarray)
    - evec[#.base, #.state] (2-dim ndarray)
    
    Note1: If a_eval_only = True, return (eval, None)
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    .      Suppose that a_args_dens[3] = N_A
    """
    l_Phi   = set_Phi_mat(a_range)
    l_Ham   = fproc_wrapper_make_Ham_mat((0, a_range, a_params, a_mass, a_args_dens, a_func_name, False))
    l_Ham_C = add_Coulomb_fpot_WS(a_range, l_Ham, a_Np, a_args_dens)
    l_ret   = eigh(l_Ham_C, l_Phi, eigvals_only=a_eval_only, eigvals=(0, a_Nstat-1))
    
    if (a_eval_only):
        return (l_ret, None)
    else:
        #for istat in range(a_Nstat): ### Normalization
        #    l_ret[1][:,istat] = l_ret[1][:,istat] / norm(l_ret[1][:,istat])
        return (l_ret[0], l_ret[1])

def make_Ham_mat_fpot_WS(a_range, a_params, a_mass, a_args_dens, a_func_name, a_Nproc = 1):
    """
    The function to make the Hamiltonian matrixes for Gauss expansion method.
    
    For arguments,
    - a_range [#.base]          (1-dim ndarray)
    - a_params[#.conf, #.param] (2-dim ndarray)
    - a_args_dens[4]            (1-dim ndarray)
    
    return: Ham[#.conf, #.base, #.base] (3-dim ndarray)
    
    Note: Suppose that a_args_dens[0] = rho_0
    .     Suppose that a_args_dens[1] = R_A
    .     Suppose that a_args_dens[2] = a_A
    .     Suppose that a_args_dens[3] = N_A
    """
    l_Nconf = len(a_params[:,0])
    
    print("#\n# Make Hamiltonian matrixes for GEM...")
    if (a_Nproc == 1):
        r_Ham = array([fproc_wrapper_make_Ham_mat((iconf, a_range, a_params[iconf,:], a_mass, 
                                                   a_args_dens, a_func_name, True)) for iconf in range(l_Nconf)])
    else:
        args_procs = [(iconf, a_range, a_params[iconf,:], a_mass, a_args_dens, a_func_name, True) 
                      for iconf in range(l_Nconf)]
        with Pool(a_Nproc) as proc:
            r_Ham = array(proc.map(fproc_wrapper_make_Ham_mat, args_procs))
    print("# Make Hamiltonian matrixes for GEM... all end\n#")
    
    return r_Ham

def solve_sch_GEM_fpot_WS_from_Ham_mat(a_Nstat, a_range, a_Ham_mat, a_args_dens, a_Np = 0, a_eval_only = True, a_Nproc = 1):
    """
    The function to solve the Schrodinger equation by using Gauss expansion method.
    ~~~ From Hamiltonian matrixes ~~~
    
    For arguments,
    - a_range[#.base]                 (1-dim ndarray)
    - a_Ham  [#.conf, #.base, #.base] (3-dim ndarray)
    - a_args_dens[3]                  (1-dim ndarray)
    
    return: (eval, evec)
    - eval[#.conf, #.state]         (2-dim ndarray)
    - evec[#,conf, #.base, #.state] (3-dim ndarray)
    
    Note1: If a_eval_only = True, return (eval, None)
    Note2: Suppose that a_args_dens[0] = rho_0
    .      Suppose that a_args_dens[1] = R_A
    .      Suppose that a_args_dens[2] = a_A
    .      Suppose that a_args_dens[3] = N_A
    """
    l_Nbase = len(a_range)
    l_Nconf = len(a_Ham_mat[:,0,0])
    l_Phi   = set_Phi_mat(a_range)
    l_Ham_C = array([add_Coulomb_fpot_WS(a_range, a_Ham_mat[iconf,:,:], a_Np, a_args_dens) for iconf in range(l_Nconf)])
    
    print("#\n# Solve Schrodinger equation by GEM...")
    if (a_Nproc == 1):
        l_ret = array([fproc_wrapper_solve_sch_GEM((iconf, l_Ham_C[iconf,:,:], l_Phi, a_eval_only, a_Nstat, True))
                       for iconf in range(l_Nconf)])
    else:
        args_procs = [(iconf, l_Ham_C[iconf,:,:], l_Phi, a_eval_only, a_Nstat, True) for iconf in range(l_Nconf)]
        with Pool(a_Nproc) as proc:
            l_ret = array(proc.map(fproc_wrapper_solve_sch_GEM, args_procs))
    print("# Solve Schrodinger equation by GEM... all end\n#")
    
    if (a_eval_only):
        return (l_ret, None);        
    else:
        return (array([l_ret[iconf, 0] for iconf in range(l_Nconf)]),
                array([l_ret[iconf, 1] for iconf in range(l_Nconf)]))
