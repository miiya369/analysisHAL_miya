# -*- coding: utf-8 -*-

"""The module for miscellaneous functions."""

def frange(a_start, a_end, a_step):
    """The function of float version of range()"""

    n = a_start
    while (n + a_step < a_end):
        yield n
        n += a_step
