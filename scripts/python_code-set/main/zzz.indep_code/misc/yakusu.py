#!/usr/bin/python

from __future__ import print_function

if __name__ == "__main__":
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 2):
        exit("usage: python %s [value]" % basename(argv[0]))
    
    val = int(argv[1])
    for i in range(2, val):
        if (val % i == 0):
            print("%d " % i, end="")
    print()
