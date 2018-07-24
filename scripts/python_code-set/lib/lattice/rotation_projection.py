# -*- coding: utf-8 -*-

"""The module to calculate the cubic group transformation."""

from numpy import roll

def rot_proj(a_data, a_state = "A1"):
    """
    The function to calculate the cubic group transformation.
    
    For arguments,
    - a_data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    - a_state (= 'A1' or 'A2' or 'E' or 'T1' or 'T2')
    
    return: r_Data[#.Zsite, #.Ysite, #.Xsite] (3-dim ndarray)
    
    Note: Assuming #.Zsite == #.Ysite == #.Xsite
    """
    if (len(a_data[:,0,0]) != len(a_data[0,:,0]) or
        len(a_data[:,0,0]) != len(a_data[0,0,:])):
        print("\nERROR: #.site (Nx,Ny,Nz) must be the same for rot_proj.\n")
        return None
    
    if   (a_state == "A1"):
        factor = 1.0 / 24.0
        l_Rchr = [1, 1, 1, 1, 1]
    elif (a_state == "A2"):
        factor = 1.0 / 24.0
        l_Rchr = [1,-1, 1, 1,-1]
    elif (a_state ==  "E"):
        factor = 2.0 / 24.0
        l_Rchr = [2, 0, 2,-1, 0]
    elif (a_state == "T1"):
        factor = 3.0 / 24.0
        l_Rchr = [3, 1,-1, 0,-1]
    elif (a_state == "T2"):
        factor = 3.0 / 24.0
        l_Rchr = [3,-1,-1, 0, 1]
    else:
        print("\nERROR: Invalid representation type (%s)\n" % a_state)
        return None
    
    return (
        # E
        l_Rchr[0] * a_data
        +
        # 6C4
        l_Rchr[1] *(roll(a_data.transpose(1,0,2)[::-1,:   ,:   ], 1, 0) +
                    roll(a_data.transpose(1,0,2)[:   ,::-1,:   ], 1, 1) +
                    roll(a_data.transpose(2,1,0)[:   ,:   ,::-1], 1, 2) +
                    roll(a_data.transpose(2,1,0)[::-1,:   ,:   ], 1, 0) +
                    roll(a_data.transpose(0,2,1)[:   ,::-1,:   ], 1, 1) +
                    roll(a_data.transpose(0,2,1)[:   ,:   ,::-1], 1, 2))
        +
        # 3C2
        l_Rchr[2] *(roll(a_data[:   ,::-1,::-1], 1, (1,2)) +
                    roll(a_data[::-1,:   ,::-1], 1, (0,2)) +
                    roll(a_data[::-1,::-1,:   ], 1, (0,1)))
        +
        # 8C3
        l_Rchr[3] *(a_data.transpose(2,0,1) +
                    roll(a_data.transpose(2,0,1)[:   ,::-1,::-1], 1, (1,2)) +
                    roll(a_data.transpose(2,0,1)[::-1,:   ,::-1], 1, (0,2)) +
                    roll(a_data.transpose(2,0,1)[::-1,::-1,:   ], 1, (0,1)) +
                    a_data.transpose(1,2,0) +
                    roll(a_data.transpose(1,2,0)[:   ,::-1,::-1], 1, (1,2)) +
                    roll(a_data.transpose(1,2,0)[::-1,:   ,::-1], 1, (0,2)) +
                    roll(a_data.transpose(1,2,0)[::-1,::-1,:   ], 1, (0,1)))
        +
        # 6C2
        l_Rchr[4] *(roll(a_data.transpose(1,0,2)[:   ,:   ,::-1], 1,  2     ) +
                    roll(a_data.transpose(1,0,2)[::-1,::-1,::-1], 1, (0,1,2)) +
                    roll(a_data.transpose(2,1,0)[:   ,::-1,:   ], 1,  1     ) +
                    roll(a_data.transpose(2,1,0)[::-1,::-1,::-1], 1, (0,1,2)) +
                    roll(a_data.transpose(0,2,1)[::-1,:   ,:   ], 1,  0     ) +
                    roll(a_data.transpose(0,2,1)[::-1,::-1,::-1], 1, (0,1,2)))
        ) * factor
