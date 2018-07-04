#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../lib")
import numpy as np
import time

### ================== Global Parameters Init. ================= ###
ifname = None
### =========================== Main =========================== ###

def main():
    return 0

### ============================================================ ###
###### Functions for arguments
def usage(ARGV0):
    print("usage  : python %s [ifile] {options}\n" % os.path.basename(ARGV0))
    print("options:")
    print("      --none [None]")
    exit(1)

def check_args():
    print("# === Check Arguments ===")
    print("# ifile =", ifname)
    print("# =======================")

def set_args(ARGC, ARGV):
    global ifname
    
    if (ARGV[1][0] == '-'):
        usage(ARGV[0])
    
    ifname = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (len(ARGV[i]) == 1):
            continue
        if (ARGV[i][0] == '-' and ARGV[i][1] == '-'):
            if   (ARGV[i] == '--pass'):
                pass
            elif (ARGV[i] == '--pass'):
                pass
            else:
                print("\nERROR: Invalid option '%s'\n" % ARGV[i])
                usage(ARGV[0])
    
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
