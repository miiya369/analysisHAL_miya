#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname  = None
ofname  = "/dev/null"
mass    = 500.0
max_r   = 10.0
range_a = 0.8
Nbase   = 10
Nstat   = 1
Np      = 0
Nproc   = 1

do_stochastic = False
eval_only     = True
from_Ham_mat  = False

def fproc_wrapper(args):
    from sch_gauss_exp.solve_sch_GEM import solve_sch_GEM
    a_iconf      = args[0]
    a_Nstat      = args[1]
    a_ranges     = args[2]
    a_params     = args[3]
    a_mass       = args[4]
    a_func_name  = args[5]    
    a_Np         = args[6]
    a_eval_only = args[7]
    results = solve_sch_GEM(a_Nstat, a_ranges, a_params, a_mass, a_func_name, a_Np, a_eval_only)
    print("# Solve Schrodinger equation by GEM... %3d end" % a_iconf)
    return results
### =========================== Main =========================== ###

def main():
    from multiprocessing import Pool
    
    from fitting.io_params           import input_params
    from sch_gauss_exp.print_eigen   import print_eigen
    from sch_gauss_exp.solve_sch_GEM import solve_sch_GEM, make_Ham_mat, solve_sch_GEM_from_Ham_mat
    from sch_gauss_exp.io_Ham_mat    import input_Ham_mat, output_Ham_mat
    
### Input parameters ###
    if (ifname != 'DEBUG'):
        if (from_Ham_mat):
            ### Calculation from Hamiltonian matrix
            max_r_in, ranges, Ham = input_Ham_mat(ifname)
            if (ranges is Ham is None):
                return -1
            eigens = solve_sch_GEM_from_Ham_mat(Nstat, ranges, Ham, Np, eval_only, Nproc)
            print_eigen(eigens[0], eigens[1], ranges, max_r_in, 0.01, eval_only)
            return 0
        else:
            func_name, params = input_params(ifname)
            if (func_name is params is None):
                return -1
            Nconf = len(params[:,0])
    else:
        ### For Debug
        from misc_QM.scattering_Idogata import Eb_Idogata3D_Swave
        global mass
        #Vzero     = -50.0; Rzero = 2.0; mass = 500.0
        Vzero     = -100.0; Rzero = 4.0; mass = 1000.0
        func_name = "SW"
        Nconf     = 1
        params    = np.array([[Vzero, Rzero]])
        print("# DEGUB MODE...")
        print("# V_0 (MeV) =", Vzero)
        print("# R_0 ( fm) =", Rzero)
        print("# mass(MeV) =", mass)
        print("# Func type = 3-dim Idogata")
        print("# Expect En =", Eb_Idogata3D_Swave(Vzero, Rzero, mass)))
    
### Construct the ranges ###
    if (not do_stochastic):
        ranges = np.array([max_r * pow(range_a, i) for i in range(Nbase)])
    else:
        ### The stocastic variational method
        ranges = np.empty(0); tmp_eval = 1e+30
        for i in range(Nbase):
            ranges = np.append(ranges, 0.0)
            while(True):
                ranges[i] = np.random.random() * max_r
                eval = solve_sch_GEM(1, ranges, params[0,:], mass, func_name, Np, True)[0]
                if (eval < tmp_eval):
                    tmp_eval = eval
                    break
    
### Solve the Schrodinger equation by GEM ###
    if (False): ### Two ways to solve GEM (first one will be deleted in future)
        print("#\n# Solve Schrodinger equation by GEM...")
        if (Nproc == 1):
            eigens = np.array([fproc_wrapper((iconf, Nstat, ranges, params[iconf,:], mass, func_name, Np, eval_only))
                               for iconf in range(Nconf)])
        else:
            args_procs = [(iconf, Nstat, ranges, params[iconf,:], mass, func_name, Np, eval_only) for iconf in range(Nconf)]
            with Pool(Nproc) as proc:
                eigens = np.array(proc.map(fproc_wrapper, args_procs))
        
        eigens = (np.array([eigens[iconf,0] for iconf in range(Nconf)]),
                  np.array([eigens[iconf,1] for iconf in range(Nconf)]))
        print("# Solve Schrodinger equation by GEM... all end\n#")
    else:
        Ham = make_Ham_mat(ranges, params, mass, func_name, Nproc)
        output_Ham_mat(ofname, max_r, ranges, Ham)
        
        eigens = solve_sch_GEM_from_Ham_mat(Nstat, ranges, Ham, Np, eval_only, Nproc)
    
### print results ###
    print_eigen(eigens[0], eigens[1], ranges, max_r, 0.01, eval_only)
    
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile (fit-param data)] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --ofile       [Output file name for matrix] Default =", ofname)
    print("      --mass        [reduced mass          (MeV)] Default =", mass)
    print("      --Np          [#.charge                   ] Default =", Np)
    print("      --max_r       [Max range of gaussian ( fm)] Default =", max_r)
    print("      --range_a     [The base of range          ] Default =", range_a)
    print("      --Nbase       [#.base                     ] Default =", Nbase)
    print("      --Nstat       [#.eigen state to calculate ] Default =", Nstat)
    print("      --Nproc       [#.process                  ] Default =", Nproc)
    print("      --stochastic   Using stochastic variational method")
    print("      --out_wave     Output wave function")
    print("      --from_mat     Input Hamiltonian matrixes")
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile     =", ifname)
    print("# ofile     =", ofname)
    print("# mass      =", mass)
    print("# Np        =", Np)
    print("# max range =", max_r)
    print("# range a   =", range_a)
    print("# N.base    =", Nbase)
    print("# N.state   =", Nstat)
    print("# N.process =", Nproc)
    print("# do_stoc.  =", do_stochastic)
    print("# Out wave  =", (not eval_only))
    print("# From mat  =", from_Ham_mat)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname, ofname, mass, max_r, range_a, Nbase, Nstat, Np
    global do_stochastic, eval_only, from_Ham_mat, Nproc
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--mass'):
                mass = float(ARGV[i+1])
            elif (ARGV[i] == '--max_r'):
                max_r = float(ARGV[i+1])
            elif (ARGV[i] == '--ofile'):
                ofname = ARGV[i+1].strip()
            elif (ARGV[i] == '--range_a'):
                range_a = float(ARGV[i+1])
            elif (ARGV[i] == '--Nbase'):
                Nbase = int(ARGV[i+1])
            elif (ARGV[i] == '--Nstat'):
                Nstat = int(ARGV[i+1])
            elif (ARGV[i] == '--Np'):
                Np = int(ARGV[i+1])
            elif (ARGV[i] == '--Nproc'):
                Nproc = int(ARGV[i+1])
            elif (ARGV[i] == '--out_wave'):
                eval_only = False
            elif (ARGV[i] == '--stochastic'):
                do_stochastic = True
            elif (ARGV[i] == '--from_mat'):
                from_Ham_mat = True
            else:
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
    if (do_stochastic):
        range_a       = "stochastic"
    
    if (from_Ham_mat and ifname != 'DEBUG'):
        ofname        = "/dev/null"
        mass          = "From mat"
        max_r         = "From mat"
        range_a       = "From mat"
        Nbase         = "From mat"
        do_stochastic = "From mat"
    
    check_args()
### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    if (argc == 1):
        usage(argv[0])
    
    set_args(argc, argv)
    
    t_start = time.time()
    if (main() != 0):
        exit("ERROR EXIT.")
    print("#\n# Elapsed time [s] = %d" % (time.time() - t_start))
