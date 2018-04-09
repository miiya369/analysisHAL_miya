#!/usr/bin/python
# -*- coding: utf-8 -*-

def frange(a_start, a_end, a_step):
    """The function of float version of range()"""

    n = a_start
    while (n + a_step < a_end):
        yield n
        n += a_step

if (__name__ == "__main__"):
    import sys
    argv = sys.argv; argc = len(argv)
    
    if (argc != 5):
        print("usage: %s [reduced mass] [Pmin] [Pdel] [Pmax]" % argv[0])
        quit()
    
    Mu   = float(argv[1].strip())
    Pmin = float(argv[2].strip())
    Pdel = float(argv[3].strip())
    Pmax = float(argv[4].strip())
    
    for p in frange(Pmin, Pmax, Pdel):
        print p**2 / (2.0 * Mu)
