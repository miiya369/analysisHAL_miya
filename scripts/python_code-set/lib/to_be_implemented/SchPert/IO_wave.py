# -*- coding: utf-8 -*-

"""The module to input the parameters of wave function."""

from numpy import array

def input_wave(a_iFname):
    """
    The function to input the parameters of wave function.
    
    return: (Range[#.base], Coeff[#.base, #.conf])
    """
    
    for Line in open(a_iFname, 'r'):
        line = Line.split()
        
        if (len(line) < 3):
            continue
        
        elif (line[1] == 'N.base'):
            if (line[3].strip().isdigit()):
                l_Nbase = int(line[3].strip())
                r_Range = array([None]*l_Nbase)
        
        elif (line[1] == 'N.conf'):
            if (line[3].strip().isdigit()):
                l_Nconf = int(line[3].strip())
                r_Coeff = array([[None]*l_Nconf]*l_Nbase)
        
        elif (line[1] == 'Range'):
            if(l_Nbase != len(line) - 3):
                print("\nERROR: #.base is differ, %d != %d, exit.\n" % (l_Nbase, len(line) - 3))
                return (None, None)
            
            for i in range(l_Nbase):
                r_Range[i] = float(line[3 + i].strip())
        
        elif (line[1] == 'Coeff'):
            if(l_Nbase != len(line) - 4):
                print("\nERROR: #.base is differ, %d != %d, exit.\n" % (l_Nbase, len(line) - 4))
                return (None, None)
            
            if (line[2].strip().isdigit()):
                iconf = int(line[2].strip())
            
            for i in range(l_Nbase):
                r_Coeff[i, iconf] = float(line[4 + i].strip())
        
        elif (Line == '# Return 1'):
            break
    
    for ibase in range(l_Nbase):
        if (r_Range[ibase] is None):
            print("\nERROR: Input error, exit.\n")
            return (None, None)
        for iconf in range(l_Nconf):
            if (r_Coeff[ibase, iconf] is None):
                print("\nERROR: Input error, exit.\n")
                return (None, None)
    
    print("# Successful to input the parameters of wave functions")
    print("# N.base  ="),; print l_Nbase
    print("# N.conf  ="),; print l_Nconf
    
    return (r_Range, r_Coeff)
