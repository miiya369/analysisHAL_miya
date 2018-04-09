#!/usr/bin/python

Nconf  = 1
Tmin   = 0
Tmax   = 1
Trange = 1

def SetPath(CONF, HAD, TSHIFT, XYZSHIFT = "000.000.000"): # <- You may change XYZshift
    return "%s/%s_correlator.+%s.%s.%s" % (CONF, HAD, TSHIFT, XYZSHIFT, CONF),

def ReadCorr(iBASE, PS_SS, CONFS, Hadron, TSHIFT, CORR):
    print
    print(" Read: correlator.%s.dir %s START" % (PS_SS, Hadron))
    for i in range(Nconf):
        ifname = "%s/correlator.%s.dir/%s" % (iBASE, PS_SS, SetPath(CONFS[i].strip(), Hadron, Tshift)[0])
        fcorr  = open(ifname, 'r')
        Lcorr  = fcorr.readlines()
        fcorr.close()
        for t in range(Trange):
            Wcorr = Lcorr[t+Tmin].split()
            CORR[i + Nconf * t] = float(Wcorr[1])
    print(" Read: correlator.%s.dir %s END" % (PS_SS, Hadron))

def MakeJK(CorrOrg, CorrJK, do_JK = 0): # <- You may change do_JK (1 = Do Jack Knife mean)
    if (Nconf != 1 and do_JK == 1):
        for t in range(Trange):
            TmpSum = 0.0
            
            for i in range(Nconf):
                TmpSum += CorrOrg[i + Nconf * t]
            for i in range(Nconf):
                CorrJK[i + Nconf * t] = (TmpSum - CorrOrg[i + Nconf * t]) / float(Nconf-1)
    else:
        for l in range(Trange * Nconf):
            CorrJK[l] = Corr[l]

def MakeErrors(DATA, ERROR, is_JKdata = 1): # <- You may change is_JKdata (1 = make Jack Knife error) 
    if (is_JKdata == 1):
        factor = float(Nconf-1)
    else:
        factor = 1.0 / float(Nconf-1)
    
    MEAN  = [0.0] * Trange
    MEAN2 = [0.0] * Trange
    for t in range(Trange):
        for i in range(Nconf):
            MEAN [t] += DATA[i + Nconf * t]
            MEAN2[t] += DATA[i + Nconf * t] ** 2
        
        MEAN [t] /= float(Nconf)
        MEAN2[t] /= float(Nconf)
        ERROR[t]  = math.sqrt(factor * (MEAN2[t] - MEAN[t]**2))

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
    
    print("%f +/- %f" % (MEAN, ERROR)),

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
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 10):
        print("usage: %s [conf list] [results dir] [Hadron1] [mass1] [Hadron2] [mass2] [Tshift] [t min] [t max]" % argv[0])
        print("     :")
        print("     :    mass should be Lattice Unit.")
        quit()
    
    ConfList = argv[1]
    ibase    = argv[2]
    Hadron1  = argv[3]
    mass1    = float(argv[4])
    Hadron2  = argv[5]
    mass2    = float(argv[6])
    Tshift   = argv[7]
    Tmin     = int(argv[8])
    Tmax     = int(argv[9])
    
    fconf = open(ConfList.strip(), 'r')
    Confs = fconf.readlines()
    Nconf = len(Confs)
    fconf.close()
    
    print(" Results dir    = %s" % ibase)
    print(" Conf list      = %s" % ConfList)
    print(" #.conf         = %d" % Nconf)
    print(" Time shift     = %s" % Tshift)
    print(" Hadron Name 1  = %s" % Hadron1)
    print(" Mass 1         = %f" % mass1)
    print(" Hadron Name 2  = %s" % Hadron2)
    print(" Mass 2         = %f" % mass2)
    print(" T range        = %d - %d" % (Tmin, Tmax))
    
    Trange = Tmax - Tmin + 1
    Corr   = [0.0] * Trange * Nconf
    JKcorr = [0.0] * Trange * Nconf
    Err    = [0.0] * Trange
    paramA = [0.0] * Nconf
    paramB = [0.0] * Nconf
    chisqA = [0.0] * Nconf
    chisqB = [0.0] * Nconf
    
    Zfactor1 = [0.0] * Nconf
    Zfactor2 = [0.0] * Nconf
    
    ReadCorr  (ibase, "PS", Confs, Hadron1, Tshift, Corr)
    MakeJK    (Corr  , JKcorr)
    MakeErrors(JKcorr, Err   )    
    FitExp    (JKcorr, mass1, Err, paramA, chisqA)
    
    ReadCorr  (ibase, "SS", Confs, Hadron1, Tshift, Corr)
    MakeJK    (Corr  , JKcorr)
    MakeErrors(JKcorr, Err   )
    FitExp    (JKcorr, mass1, Err, paramB, chisqB)
    
    print
    print(" Chisq/dof for PS (%s) = " % Hadron1),
    PrintMeanError(chisqA)
    print
    print(" Chisq/dof for SS (%s) = " % Hadron1),
    PrintMeanError(chisqB)
    print

    for i in range(Nconf):
        Zfactor1[i] = math.sqrt(2.0*float(mass1)) * paramA[i] / math.sqrt(paramB[i])
    
    ReadCorr  (ibase, "PS", Confs, Hadron2, Tshift, Corr)
    MakeJK    (Corr  , JKcorr)
    MakeErrors(JKcorr, Err   )
    FitExp    (JKcorr, mass2, Err, paramA, chisqA)

    ReadCorr  (ibase, "SS", Confs, Hadron2, Tshift, Corr)
    MakeJK    (Corr  , JKcorr)
    MakeErrors(JKcorr, Err   )
    FitExp    (JKcorr, mass2, Err, paramB, chisqB)

    print
    print(" Chisq/dof for PS (%s) = " % Hadron2),
    PrintMeanError(chisqA)
    print
    print(" Chisq/dof for SS (%s) = " % Hadron2),
    PrintMeanError(chisqB)
    print

    for i in range(Nconf):
        Zfactor2[i] = math.sqrt(2.0*float(mass2)) * paramA[i] / math.sqrt(paramB[i])
    
    Zfactor12 = [0.0] * Nconf
    for i in range(Nconf):
        Zfactor12[i] = Zfactor1[i] / Zfactor2[i]
    
    print
    print(" sqrt(Z_1/Z_2) (Point)  = "),
    PrintMeanError(Zfactor12)
