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

def PrintFunc(FUNC_TYPE, PARAM): # <- CHANGE HERE !
    if (FUNC_TYPE == 2):
        print("# %f*exp(-(%f)*x)"       % (PARAM[0], PARAM[1]))
    elif (FUNC_TYPE == 3):
        print("# %f*exp(-(x/(%f))**2)"  % (PARAM[0], PARAM[1]))
    elif (FUNC_TYPE == 4):
        print("# %f*exp(-(x/(%f))**2)+" % (PARAM[0], PARAM[1]) +
              "(%f)*exp(-(x/(%f))**2)"  % (PARAM[2], PARAM[3]))
    elif (FUNC_TYPE == 5):
        print("# %f*exp(-(x/(%f))**2)+" % (PARAM[0], PARAM[1]) +
              "(%f)*exp(-(x/(%f))**2)+" % (PARAM[2], PARAM[3]) +
              "(%f)*exp(-(x/(%f))**2)"  % (PARAM[4], PARAM[5]))
    elif (FUNC_TYPE == 8):
        print("# %f*exp(-(x/(%f))**2)+" % (PARAM[0], PARAM[1]) +
              "(%f)*exp(-(x/(%f))**2)+" % (PARAM[2], PARAM[3]) +
              "(%f)*exp(-(x/(%f))**2)+" % (PARAM[4], PARAM[5]) +
              "(%f)*exp(-(x/(%f))**2)"  % (PARAM[6], PARAM[7]))
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
    print("# Func type = %s\n#" % Fname)
    
    params = [0.0] * (Nconf * Nparam)

    for i in range(Nconf * Nparam):
        params[i] = struct.unpack('<d', ifile.read(8))[0]
    
    ifile.close()
    
    factor = Nconf - 1
    Mean   = [0.0] * Nparam
    for i in range(Nparam):
        Mean2 = 0.0
        for j in range(Nconf):
            Mean[i] += params[j + Nconf*i]
            Mean2   += params[j + Nconf*i] ** 2
        Mean[i] /= float(Nconf)
        Mean2   /= float(Nconf)
        Error    = math.sqrt(abs(factor * (Mean2 - Mean[i] ** 2)))
        print("# param[%d] = %f +/- %f" % (i, Mean[i], Error))

    print("#\n#"),
    for i in range(Nparam):
        print("%f" % Mean[i]),
    
    print
    print("#")
    PrintFunc(funcN, Mean)
    print("#\n# r | Func Mean | Func Error | Func Min | Func Max")
    
    Npoint     = int(max_r / div_r)
    tmp_params = [0.0] * Nparam
    for i in range(Npoint):
        r = i * div_r
        
        tmpV_max = -1e10
        tmpV_min =  1e10
        
        Mean  = 0.0
        Mean2 = 0.0
        
        for j in range(Nconf):
            for k in range(Nparam):
                tmp_params[k] = params[j + Nconf*k]
            tmpV = FitFunc(r/lat_space, funcN, tmp_params) * hc_api
            
            if (tmpV_max < tmpV):
                tmpV_max = tmpV
            if (tmpV_min > tmpV):
                tmpV_min = tmpV
            
            Mean  += tmpV
            Mean2 += tmpV ** 2
        
        Mean  /= float(Nconf)
        Mean2 /= float(Nconf)
        Error  = math.sqrt(factor * (Mean2 - Mean ** 2))
        print("%f %f %f %f %f" % (r, Mean, Error, tmpV_min, tmpV_max))
