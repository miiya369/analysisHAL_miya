#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
import numpy as np

from lattice.rotation_projection import rot_proj
from struct import unpack

if (__name__ == "__main__"):
    eig_val = np.empty( 10 )
    eig_vec = np.empty((10,32,32,32))
    with open("/Users/miiya/Dropbox/programs/analysisHAL_miya/scripts"+
              "/python_code-set/zzz.main_sub/data_test/solve_sch/eig_N10.out", "rb") as ifile:
        for i in range(10):
            eig_val[i] = unpack('>d', ifile.read(8))[0]
            for z in range(32):
                for y in range(32):
                    for x in range(32):
                        eig_vec[i,z,y,x] = unpack('>d', ifile.read(8))[0]

    for z in range(32//2):
        for y in range(32//2):
            for x in range(32//2):
                r = np.sqrt(x**2+y**2+z**2) * 0.1
                print(r, end=" ")
                for i in range(10):
                    print(eig_vec[i,z,y,x], end=" ")
                print()
    print("\n")
    print("#Rotating...", end = " ")
    eig_vec = np.array([rot_proj(eig_vec[i,:,:,:]) for i in range(10)])
    print("#DONE")
    
    for z in range(32//2):
        for y in range(32//2):
            for x in range(32//2):
                r = np.sqrt(x**2+y**2+z**2) * 0.1
                print(r, end=" ")
                for i in range(10):
                    print(eig_vec[i,z,y,x], end=" ")
                print()
    print("#\n#")
    for i in range(10):
        print("#", eig_val[i] * 197.327053 / 0.1)        
