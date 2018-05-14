#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 3):
        exit("usage: python %s [value] [#.yakusu]" % basename(argv[0]))
    
    valu    = int(argv[1])
    Nyakusu = int(argv[2])
    tmp_j   = 2
    
    for i in range(Nyakusu+1):
        for j in range(tmp_j, valu):
            if (valu%j == 0):
                print("%d " % j, end="")
                tmp_j = j+1; break
    print()
