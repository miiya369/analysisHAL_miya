#!/usr/bin/python

Nparam = 1
Nconf  = 1
Fname  = "UNKNOWN"

def SetFunc(FUNC_TYPE): # <- CHANGE HERE !
    global Nparam, Fname

    if (FUNC_TYPE == 2):
        Nparam = 2
        Fname  = "one Exponential"
    elif (FUNC_TYPE == 3):
        Nparam = 2
        Fname  = "one Gaussian"
    elif (FUNC_TYPE == 4):
        Nparam = 4
        Fname  = "2-Ranges Gaussian"
    elif (FUNC_TYPE == 5):
        Nparam = 6
        Fname  = "3-Ranges Gaussian"
    elif (FUNC_TYPE == 8):
        Nparam = 8
        Fname  = "4-Ranges Gaussian"
    else:
        print("Invalid fit function type: %d" % FUNC_TYPE)
        quit()

def FitFunc(R, FUNC_TYPE, PARAM): # <- CHANGE HERE !
    if (FUNC_TYPE == 2):
        ret = 0.0
        for i in range(0, Nparam, 2):
            EX   = math.exp(-PARAM[i+1]*R)
            ret += RAPAM[i]*EX
        return ret
    elif (FUNC_TYPE == 3 or FUNC_TYPE == 4 or FUNC_TYPE == 5 or FUNC_TYPE == 8):
        ret = 0.0
        for i in range(0, Nparam, 2):
            ARG  = R/PARAM[i+1]
            EX   = math.exp(-(ARG*ARG))
            ret += PARAM[i]*EX
        return ret
    else:
        print("Invalid fit function type: %d" % FUNC_TYPE)
        quit()

if __name__ == "__main__":
    import sys
    import math
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 6):
        print "usage: %s [Fit param] [r min] [r max] [div r] [LatUnit=0, PhycUnit=LatSpacing]" % argv[0]
        quit()
    
    ifile = open(argv[1].strip(), 'rb')
    min_r = float(argv[2])
    max_r = float(argv[3])
    div_r = float(argv[4])
    
    if (float(argv[5]) == 0.0):
        lat_space = 1.0
        hc_api    = 1.0
    else:
        lat_space = float(argv[5])
        hc_api    = 197.327053 / lat_space
    
    Nconf = struct.unpack('<i', ifile.read(4))[0]
    funcN = struct.unpack('<i', ifile.read(4))[0]
    
    SetFunc(funcN)
    
    print("# N.conf    = %d" % Nconf)
    print("# N.param   = %d" % Nparam)
    print("# Func type = %s" % Fname)
    
    params = [0.0] * (Nconf * Nparam)

    for i in range(Nconf * Nparam):
        params[i] = struct.unpack('<d', ifile.read(8))[0]
    
    ifile.close()
    
    print("#\n# r | Conf1 | Conf2 | ...")
    
    Npoint     = int(max_r / div_r)
    tmp_params = [0.0] * Nparam
    for i in range(Npoint):
        r = i * div_r
        print("%f" % r),
        
        for j in range(Nconf):
            for k in range(Nparam):
                tmp_params[k] = params[j + Nconf*k]
            tmpV = FitFunc(r/lat_space, funcN, tmp_params) * hc_api
            print("%f" % tmpV),
        print
