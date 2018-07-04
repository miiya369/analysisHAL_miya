# -*- coding: utf-8 -*-

"""The module to calculate NLO potential by using SVD."""

def calc_potential_NLO(a_Pot, a_Rcorr, a_LapRcorr):
    """
    The function to calculate the NLO potential by using SVD.
    
    For arguments,
    - a_Pot     [#.channel, #.channel, #.data, #.eq] (4-dim array)
    - a_Rcorr   [#.channel, #.channel, #.data, #.eq] (4-dim array)
    - a_LapRcorr[#.channel, #.channel, #.data, #.eq] (4-dim array)
    
    return: (r_Pot__LO[#.channel, #.channel, #.data], 
             r_Pot_NLO[#.channel, #.channel, #.data])
    
    Note: #.channel, #.data, #.eq are got from a_Pot
    """
    
    from numpy        import empty, dot
    from scipy.linalg import inv, svd, diagsvd
    
    l_Ndata = len(a_Pot[0, 0, :, 0])
    l_Nch   = len(a_Pot[:, 0, 0, 0])
    l_Neq   = len(a_Pot[0, 0, 0, :])
    
    if (l_Neq < 2):
        print("\nERROR: SVD cannot solve the NLO potential in #.eq < 2, exit.\n")
        return (None, None)
    
    r_Pot__LO = empty((l_Nch, l_Nch, l_Ndata))
    r_Pot_NLO = empty((l_Nch, l_Nch, l_Ndata))
    
    l_Mat   = empty((2*l_Nch, l_Nch*l_Neq))
    l_TmpP  = empty((  l_Nch, l_Nch*l_Neq))
    l_TmpR  = empty((  l_Nch, l_Nch      ))
    l_TmpLR = empty((  l_Nch, l_Nch      ))
    
    for idata in range(l_Ndata):
        for ieq in range(l_Neq):
            for ich in range(l_Nch):
                for jch in range(l_Nch):
                    l_TmpP [ich, jch + l_Nch * ieq] = a_Pot[ich, jch, idata, ieq]
                    
                    l_TmpR [ich, jch] = a_Rcorr   [ich, jch, idata, ieq]
                    l_TmpLR[ich, jch] = a_LapRcorr[ich, jch, idata, ieq]
            
            l_TmpLR_Rinv = dot(l_TmpLR, inv(l_TmpR))
            
            for ich in range(l_Nch):
                for jch in range(l_Nch):
                    l_Mat[ich        , jch + l_Nch * ieq] = 1.0 if ich == jch else 0.0
                    l_Mat[ich + l_Nch, jch + l_Nch * ieq] = l_TmpLR_Rinv[ich, jch]
        
        l_U, l_S, l_V = svd(l_Mat)
        
        l_Res = dot(dot(dot(l_TmpP, l_V.T), diagsvd(1.0/l_S, l_Nch*l_Neq, 2*l_Nch)), l_U.T)
        
        for ich in range(l_Nch):
            for jch in range(l_Nch):
                r_Pot__LO[ich, jch, idata] = l_Res[ich, jch]
                r_Pot_NLO[ich, jch, idata] = l_Res[ich, jch + l_Nch]
        
    return (r_Pot__LO, r_Pot_NLO)
