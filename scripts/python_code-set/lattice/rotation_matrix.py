# -*- coding: utf-8 -*-

"""The module for definition of rotation matrix."""

from numpy import array

def irep(a_i):
    """
    The numbaring for the rotation matrix (a_i = 0 ~ 23)
    
    This function is used i.e. 
    'rot_mat(Lsize, *irep(i))' or 'rot_chr("A1", irep(i)[0])'
    """
    
    if   (      a_i == 0):
        return ( "E" , a_i- 0)
    elif ( 1 <= a_i <  7):
        return ("6C4", a_i- 1)
    elif ( 7 <= a_i < 10):
        return ("3C2", a_i- 7)
    elif (10 <= a_i < 18):
        return ("8C3", a_i-10)
    elif (18 <= a_i < 24):
        return ("6C2", a_i-18)
    else:
        print("\nERROR: Unkown index, '%d'." % a_i)
        return (None, None)

def rot_chr(a_state, a_rep):
    if   (a_state == "A1"):
        a_istate = 0
    elif (a_state == "A2"):
        a_istate = 1
    elif (a_state ==  "E"):
        a_istate = 2
    elif (a_state == "T1"):
        a_istate = 3
    elif (a_state == "T2"):
        a_istate = 4
    else:
        print("\nERROR: Unkown state, '%s'." % a_state)
        a_istate = None
    
    if   (a_rep ==  "E" ):
        a_irep = 0
    elif (a_rep == "6C4"):
        a_irep = 1
    elif (a_rep == "3C2"):
        a_irep = 2
    elif (a_rep == "8C3"):
        a_irep = 3
    elif (a_rep == "6C2"):
        a_irep = 4
    else:
        print("\nERROR: Unkown representation, '%s'." % a_rep)
        a_irep = None
    
    # Order: E, 6C4, 3C2, 8C3, 6C2
    r_mat = [[ 1, 1, 1, 1, 1], # A1
             [ 1,-1, 1, 1,-1], # A2
             [ 2, 0, 2,-1, 0], # E
             [ 3, 1,-1, 0,-1], # T1
             [ 3,-1,-1, 0, 1]] # T2
    
    return r_mat[a_istate][a_irep]

def rot_mat(a_L, a_rep, a_rep_i):
    """
    The definition of the rotation matrix
    
    return: rot_mat[4, 4] (2-dim ndarray)
    """
    
    if   (a_rep ==  "E"  and a_rep_i == 0):
        r_mat = [[ 1, 0, 0, 0],
                 [ 0, 1, 0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C4" and a_rep_i == 0):
        r_mat = [[ 1, 0, 0, 0],
                 [ 0, 0,-1, a_L],
                 [ 0, 1, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C4" and a_rep_i == 1):
        r_mat = [[ 0, 0, 1, 0],
                 [ 0, 1, 0, 0],
                 [-1, 0, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C4" and a_rep_i == 2):
        r_mat = [[ 0,-1, 0, a_L],
                 [ 1, 0, 0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C4" and a_rep_i == 3):
        r_mat = [[ 1, 0, 0, 0],
                 [ 0, 0, 1, 0],
                 [ 0,-1, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C4" and a_rep_i == 4):
        r_mat = [[ 0, 0,-1, a_L],
                 [ 0, 1, 0, 0],
                 [ 1, 0, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C4" and a_rep_i == 5):
        r_mat = [[ 0, 1, 0, 0],
                 [-1, 0, 0, a_L],
                 [ 0, 0, 1, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "3C2" and a_rep_i == 0):
        r_mat = [[ 1, 0, 0, 0],
                 [ 0,-1, 0, a_L],
                 [ 0, 0,-1, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "3C2" and a_rep_i == 1):
        r_mat = [[-1, 0, 0, a_L],
                 [ 0, 1, 0, 0],
                 [ 0, 0,-1, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "3C2" and a_rep_i == 2):
        r_mat = [[-1, 0, 0, a_L],
                 [ 0,-1, 0, a_L],
                 [ 0, 0, 1, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 0):
        r_mat = [[ 0, 0, 1, 0],
                 [ 1, 0, 0, 0],
                 [ 0, 1, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 1):
        r_mat = [[ 0,-1, 0, a_L],
                 [ 0, 0,-1, a_L],
                 [ 1, 0, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 2):
        r_mat = [[ 0, 0,-1, a_L],
                 [ 1, 0, 0, 0],
                 [ 0,-1, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 3):
        r_mat = [[ 0,-1, 0, a_L],
                 [ 0, 0, 1, 0],
                 [-1, 0, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 4):
        r_mat = [[ 0, 1, 0, 0],
                 [ 0, 0, 1, 0],
                 [ 1, 0, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 5):
        r_mat = [[ 0, 0, 1, 0],
                 [-1, 0, 0, a_L],
                 [ 0,-1, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 6):
        r_mat = [[ 0, 1, 0, 0],
                 [ 0, 0,-1, a_L],
                 [-1, 0, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "8C3" and a_rep_i == 7):
        r_mat = [[ 0, 0,-1, a_L],
                 [-1, 0, 0, a_L],
                 [ 0, 1, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C2" and a_rep_i == 0):
        r_mat = [[ 0, 1, 0, 0],
                 [ 1, 0, 0, 0],
                 [ 0, 0,-1, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C2" and a_rep_i == 1):
        r_mat = [[ 0,-1, 0, a_L],
                 [-1, 0, 0, a_L],
                 [ 0, 0,-1, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C2" and a_rep_i == 2):
        r_mat = [[-1, 0, 0, a_L],
                 [ 0, 0,-1, a_L],
                 [ 0,-1, 0, a_L],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C2" and a_rep_i == 3):
        r_mat = [[-1, 0, 0, a_L],
                 [ 0, 0, 1, 0],
                 [ 0, 1, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C2" and a_rep_i == 4):
        r_mat = [[ 0, 0, 1, 0],
                 [ 0,-1, 0, a_L],
                 [ 1, 0, 0, 0],
                 [ 0, 0, 0, 1]]
    elif (a_rep == "6C2" and a_rep_i == 5):
        r_mat = [[ 0, 0,-1, a_L],
                 [ 0,-1, 0, a_L],
                 [-1, 0, 0, a_L],
                 [ 0, 0, 0, 1]]
    else:
        print("\nERROR: Unknown representation, rep=%s (i=%d)\n" % (a_rep, a_rep_i))
        return None
    
    return array(r_mat)
