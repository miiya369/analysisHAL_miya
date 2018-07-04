#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import FoldingPotential.FourierBesselCoeff as fb
import FoldingPotential.DensityFunc as df
import FoldingPotential.PotentialFunc as pf
import MiscFuncs.Misc as mm

import numpy as np

A  = np.array([12, 28, 40, 58, 90, 208])
Np = np.array([ 6, 14, 20, 28, 40,  82])
RA = np.array([2.1545, 3.15, 3.6, 4.094, 4.9, 6.8])
aA = np.array([0.425, 0.475, 0.523, 0.54, 0.515, 0.515])
R0 = np.array([float(df.calc_rho0(A[i], df.dens_woods_saxon, [RA[i], aA[i]])) for i in range(6)])
#print R0; quit()

for r in mm.frange(0, 20, 0.01):
    print("%lf" % r),
    for i in range(len(A)):
        print("%lf" % (df.dens_FBcoe(r, fb.params_FBcoe(A[i])))),
    print

print; print

for r in mm.frange(0, 20, 0.01):
    print("%lf" % r),
    for i in range(len(A)):
        print("%lf" % (R0[i] * Np[i]/A[i] * df.dens_woods_saxon(r, RA[i], aA[i]))),
    print

print; print

for r in mm.frange(0, 1000, 0.1):
    print("%lf" % r),
    for i in range(len(A)):
        print("%lf" % (pf.fpot_Coulomb_dens_FBcoe(r, fb.params_FBcoe(A[i])))),
    print

print; print

for r in mm.frange(0.1, 1000, 0.1):
    print("%lf" % r),
    for i in range(len(A)):
        print("%lf" % (pf.fpot_Coulomb_dens_WS(r, R0[i], RA[i], aA[i], Np[i], A[i]))),
    print
