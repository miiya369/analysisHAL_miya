#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Takaya Miyamoto
# Brief : Repartition of 'field' files
# Date  : Thu May 17 16:19:40 JST 2018

from __future__ import print_function
from   struct   import pack

import sys, os
import numpy as np

### ================== Global Parameters Init. ================= ###
file_name_base = "%s.%02d.%02d.%02d.%02d"

file_name = "psi_vec_noise"

idir = None
odir = None

Nint = 2*3*4 # (real/imag) x (#.color) x (#.dirac)
Next = 64    # external index (e.g. #.noise or #.eigen value)

Lsize = np.array([16,16,16,32])
idiv  = np.array([2,2,2,4])
odiv  = np.array([1,1,1,1])

data = None
### =========================== Main =========================== ###
def main():
    print("#\n### Main part START ###")
    
    iNsize = np.array([Lsize[i]//idiv[i] for i in range(4)])
    oNsize = np.array([Lsize[i]//odiv[i] for i in range(4)])
    
    for it in range(idiv[3]):
        for iz in range(idiv[2]):
            for iy in range(idiv[1]):
                for ix in range(idiv[0]):
                    ifname = ("%s/%s" % (idir, (file_name_base % (file_name,ix,iy,iz,it))))
                    print("# Read from", ifname)
                    
                    idata = np.fromfile(ifname, '>d').reshape(Next,iNsize[3],iNsize[2],
                                                              iNsize[1],iNsize[0],Nint)
                    data[:,
                        it*iNsize[3]:(it+1)*iNsize[3],
                        iz*iNsize[2]:(iz+1)*iNsize[2],
                        iy*iNsize[1]:(iy+1)*iNsize[1],
                        ix*iNsize[0]:(ix+1)*iNsize[0],:] = idata
    
    if (odir is None):
        return 0
    print("#");
    
    for it in range(odiv[3]):
        for iz in range(odiv[2]):
            for iy in range(odiv[1]):
                for ix in range(odiv[0]):
                    ofname = ("%s/%s" % (odir, (file_name_base % (file_name,ix,iy,iz,it))))
                    print("# Write to ", ofname)
                    
                    odata = data[:,
                        it*oNsize[3]:(it+1)*oNsize[3],
                        iz*oNsize[2]:(iz+1)*oNsize[2],
                        iy*oNsize[1]:(iy+1)*oNsize[1],
                        ix*oNsize[0]:(ix+1)*oNsize[0],:].flatten()
                    with open(ofname, 'wb') as ofdata:
                        ofdata.write(pack('>%dd' % len(odata), *odata))
    
    print("### Main part  END  ###")
    return 0

### ============================================================ ###
###### Functions for arguments ######
def usage(ARGV0):
    print("usage  : python %s [input directory] {options}" % os.path.basename(ARGV0))
    print("options:")
    print("      --fname [base of file name (@@@.ix.iy.iz.it)] Default =", file_name)
    print("      --odir  [output directory                   ] Default =", odir)
    print("      --Next  [#.external index                   ] Default =", Next)
    print("      --Nint  [#.internal index                   ] Default =", Nint)
    print("      --Lsize [global lattice size                ] Default =", Lsize)
    print("      --idiv  [#.node division for input (x,y,z,t)] Default =", idiv)
    print("      --odiv  [#.node division for input (x,y,z,t)] Default =", odiv)
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# fname =", file_name)
    print("# idir  =", idir)
    print("# odir  =", odir)
    print("# N.ext =", Next)
    print("# N.int =", Nint)
    print("# Lsize =", Lsize)
    print("# idiv  =", idiv)
    print("# odiv  =", odiv)
    print("# =======================")

def set_args(ARGC, ARGV):
    global file_name, idir, odir, Next, Nint, Lsize, idiv, odiv, data
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    idir = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--odir'):
                odir = ARGV[i+1].strip()
            elif (ARGV[i] == '--fname'):
                file_name = ARGV[i+1].strip()
            elif (ARGV[i] == '--Next'):
                Next = int(ARGV[i+1].strip())
            elif (ARGV[i] == '--Nint'):
                Nint = int(ARGV[i+1].strip())
            elif (ARGV[i] == '--Lsize'):
                Lsize[:] = np.array([int(ARGV[i+1]),
                                     int(ARGV[i+2]),
                                     int(ARGV[i+3]),
                                     int(ARGV[i+4])])
            elif (ARGV[i] == '--idiv'):
                idiv[:] = np.array([int(ARGV[i+1]),
                                    int(ARGV[i+2]),
                                    int(ARGV[i+3]),
                                    int(ARGV[i+4])])
            elif (ARGV[i] == '--odiv'):
                odiv[:] = np.array([int(ARGV[i+1]),
                                    int(ARGV[i+2]),
                                    int(ARGV[i+3]),
                                    int(ARGV[i+4])])
            else:
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
    check_args()
    
    ### Check
    for i in range(4):
        if (Lsize[i] % idiv[i] != 0):
            exit("\nERROR: Invalid size for (Lsize[%d], idiv[%d]) = (%d, %d), exit.\n" %
                 (i, i, Lsize[i], idiv[i]))
        if (Lsize[i] % odiv[i] != 0):
            exit("\nERROR: Invalid size for (Lsize[%d], odiv[%d]) = (%d, %d), exit.\n" %
                 (i, i, Lsize[i], odiv[i]))
    
    data = np.empty((Next, Lsize[3], Lsize[2], Lsize[1], Lsize[0], Nint))

### ============================================================ ###
### ============================================================ ###
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    
    if (argc == 1):
        usage(argv[0])
    
    set_args(argc, argv)
    
    if (main() != 0):
        exit("ERROR EXIT.")
