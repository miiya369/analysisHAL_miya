#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../lib")
import numpy as np

from common.misc          import frange
from common.statistics    import make_mean_err
from fitting.io_params    import input_params
from fitting.fitfunc_type import set_fitfunc_from_fname

### ================== Global Parameters Init. ================= ###
Npot    = 4
rot_mat = np.array([[1, -3, -3,  9],
                    [1,  1, -3, -3],
                    [1, -3,  1, -3],
                    [1,  1,  1,  1]])

fbase = "/Users/miiya/data/Kfit"
ifile = ["Pot.rot.opr.NX_.__0.__0.t11.fitparam.2G1Ysq.rmax.3.0",
         "Pot.rot.opr.NX_.sig.__0.t11.fitparam.2G.rmax.3.0",
         "Pot.rot.opr.NX_.__0.tau.t11.fitparam.2G.rmax.3.0",
         "Pot.rot.opr.NX_.sig.tau.t11.fitparam.2G1Y.rmax.3.0"]

r_min  = 0.001
r_del  = 0.01
r_max  = 2.5
### =========================== Main =========================== ###
def main():
    ffunc = [set_fitfunc_from_fname(input_params(fbase+"/"+ifile[i], False)[0]) for i in range(Npot)]
    param = [                       input_params(fbase+"/"+ifile[i], False)[1]  for i in range(Npot)]
    
    Nconf = len(param[0][:,0])
    ### check #.conf
    for i in range(Npot):
        if (Nconf != len(param[i][:,0])):
            print("\nERROR: Nconf is differ, exit\n"); return -1
    
    for r in frange(r_min, r_max, r_del):
        ipot = np.array([[ffunc[i](r, *param[i][iconf,:]) for i in range(Npot)] for iconf in range(Nconf)])
        rpot = np.array([np.dot(rot_mat, ipot[iconf,:])          for iconf in range(Nconf)])
        
        print("%lf %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e" % 
              (r, 
               *make_mean_err(rpot[:,0]), *make_mean_err(rpot[:,1]),
               *make_mean_err(rpot[:,2]), *make_mean_err(rpot[:,3])))
    
    return 0

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    if (main() != 0):
        exit("ERROR EXIT.")
