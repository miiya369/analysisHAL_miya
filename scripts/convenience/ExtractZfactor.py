#!/usr/bin/python

Nconf  = 1
Tmin   = 0
Tmax   = 1
Trange = 1

def PrintMeanError(DATA, is_JKdata = 1): # <- You may change is_JKdata (1 = make Jack Knife error)
    if (is_JKdata == 1):
        factor = float(Nconf-1)
    else:
        factor = 1.0 / float(Nconf-1)
    
    MEAN  = 0.0
    MEAN2 = 0.0
    for i in range(Nconf):
        MEAN  += DATA[i]
        MEAN2 += DATA[i] ** 2
    
    MEAN  /= float(Nconf)
    MEAN2 /= float(Nconf)
    ERROR  = math.sqrt(factor * (MEAN2 - MEAN**2))
    
    if (MEAN > 0.0001):
        print("%f +/- %f" % (MEAN, ERROR))
    else:
        print("%e +/- %e" % (MEAN, ERROR))

def FitExp(DATA, MASS, ERR, PARAM, CHISQ):
    for i in range(Nconf):
        tmpA     = 0.0
        tmpB     = 0.0
        CHISQ[i] = 0.0
        for t in range(Trange):
            tmpA += (                  math.exp(-2 * MASS * (t+Tmin))) / (ERR[t] ** 2)
            tmpB += (DATA[i+Nconf*t] * math.exp(-    MASS * (t+Tmin))) / (ERR[t] ** 2)
        
        PARAM[i] = tmpB / tmpA
        for t in range(Trange):
            CHISQ[i] += ((PARAM[i]*math.exp(-MASS * (t+Tmin)) - DATA[i+Nconf*t])**2) / (ERR[t] ** 2)
            CHISQ[i] /= float(Trange - 1)

if __name__ == "__main__":
    import sys
    import math
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 5 and argc != 6):
        print("usage: %s [corr.fit file] [mass (Lat.Unit)] [t_min] [t_max]" % argv[0])
        quit()
    
    Tmin = int(argv[3])
    Tmax = int(argv[4])
    
    if (Tmin < 0):
        print("t_min should be greater than 0.")
        quit()
    if (Tmin >= Tmax):
        print("t_max should be greater than t_min.")
        quit()
    
    ifile = open(argv[1].strip(), 'rb')
    mass  = float(argv[2])
    
    ftype = struct.unpack('<i', ifile.read(4))[0]
    if (ftype != 1):
        print("The input file is not 'Miyamoto-format correlator-type fit file'.")
        quit()
        
    Nconf = struct.unpack('<i', ifile.read(4))[0]
    Ndata = struct.unpack('<i', ifile.read(4))[0]
    
    if (Tmax > Ndata):
        print("ERROR: t_max(%d) > Range of time(%d)" % (Tmax, Ndata))
        quit()
    
    Trange = Tmax - Tmin + 1
    Corr   = [0.0] * Trange * Nconf
    Err    = [0.0] * Trange
    param  = [0.0] * Nconf
    chisq  = [0.0] * Nconf
    
    for it in range(Ndata):
        if (Tmin <= it and it <= Tmax):
            Err[it-Tmin] = struct.unpack('<d', ifile.read(8))[0]
        else:
            dummy = struct.unpack('<d', ifile.read(8))[0]
        
        for iconf in range(Nconf):
            if (Tmin <= it and it <= Tmax):
                Corr[iconf + Nconf * (it-Tmin)] = struct.unpack('<d', ifile.read(8))[0]
            else:
                dummy = struct.unpack('<d', ifile.read(8))[0]
    
    FitExp(Corr, mass, Err, param, chisq)

    if (argc == 6):
        print("chisq:")
        PrintMeanError(chisq)
        print("param:")
        PrintMeanError(param)
    else:
        for iconf in range(Nconf):
            print("%e" % param[iconf])
