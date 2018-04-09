#!/usr/bin/python

### Auther: Takaya Miyamoto
### Date  : Fri Oct  7 20:11:23 JST 2016
### Brief : Chack the bin size dependence of Jack Knife sample

Nconf = 1
Ntime = 1

col_corr = 0 # <- The column of correlator-file to read 
             # (For complex, 0 => column 1 & 2, 1 => column 3 & 4, ... are read)

def makeJK(DATAorg, DATAjk, bin_size):
    Nbin   = Nconf / bin_size
    TmpSum = complex(0.0, 0.0)
    
    for i in range(Nconf):
        TmpSum += DATAorg[i]
    
    for i in range(Nbin):
        TmpSum2 = TmpSum
        
        for j in range(bin_size):
            TmpSum2 -= DATAorg[j + bin_size * i]
            
        DATAjk[i] = TmpSum2 / float(Nconf-bin_size)

def makeJKerr(DATA, bin_size):
    Nbin   = Nconf / bin_size
    
    MeanR  = 0.0
    Mean2R = 0.0
    MeanI  = 0.0
    Mean2I = 0.0
    
    for ibin in range(Nbin):
        MeanR  += DATA[ibin].real
        Mean2R += DATA[ibin].real ** 2
        MeanI  += DATA[ibin].imag
        Mean2I += DATA[ibin].imag ** 2
    
    MeanR  /= float(Nbin)
    Mean2R /= float(Nbin)
    MeanI  /= float(Nbin)
    Mean2I /= float(Nbin)
    ErrR    = math.sqrt((Nbin-1) * (Mean2R - MeanR ** 2))
    ErrI    = math.sqrt((Nbin-1) * (Mean2I - MeanI ** 2))
    
    print "%d %e %e %e %e" % (bin_size, MeanR, ErrR, MeanI, ErrI)

if (__name__ == "__main__"):
    import sys
    import math
    import cmath
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc < 5):
        print("usage: %s [conf-list] [corr-file] [0=OrgCorr, 1=Effmass] [t1] [t2] ..." % argv[0])
        print("     : At the corr-file, \"REPCONF\" will be replaced by conf-name in the conf-list.")
        quit()
    
    ConfList  = argv[1]
    CorrFile  = argv[2]
    CorEffFlg = int(argv[3])
    Ntime     = argc - 4
    
    if (CorEffFlg != 0 and CorEffFlg != 1):
        print("usage: %s [conf-list] [corr-file] [0=OrgCorr, 1=Effmass] [t1] [t2] ..." % argv[0])
        print("     : At the corr-file, \"REPCONF\" will be replaced by conf-name in the conf-list.")
        quit()
    
    fconf = open(ConfList.strip(), 'r')
    Confs = fconf.readlines()
    Nconf = len(Confs)
    fconf.close()
    
    if (Nconf < 2):
        print("#.conf should be > 1.")
        quit()
    
    fcorr = open(CorrFile.replace('REPCONF', Confs[0].strip()), 'r')
    Lcorr = fcorr.readlines()
    Tsize = len(Lcorr)
    fcorr.close()
    
    print("# Setup:")
    print("#    - Conf list = %s" % ConfList)
    print("#    - Corr file = %s" % CorrFile)
    print("#    - #.conf    = %d" % Nconf)
    print("#    - time size = %d" % Tsize)
    print("#    - #.time    = %d" % Ntime)
    if (CorEffFlg == 0):
        print("#    - Check for : correlator")
    else:
        print("#    - Check for : effective mass")
    print("#")
    
    Corr_p0 = [0.0] * Nconf
    
    if (CorEffFlg == 1):
        Corr_p1 = [0.0] * Nconf
    
    for it in range(Ntime):
        time = int(argv[it+4])
        print("### time slice = %d" % time)
        print
        
        for iconf in range(Nconf):
            fcorr = open(CorrFile.replace('REPCONF', Confs[iconf].strip()), 'r')
            Lcorr = fcorr.readlines()
            fcorr.close()
            
            Wcorr = Lcorr[time].split()
            Corr_p0[iconf] = complex(float(Wcorr[1+2*col_corr]), float(Wcorr[2+2*col_corr]))
            
            if (CorEffFlg == 1):
                Wcorr = Lcorr[time+1].split()
                Corr_p1[iconf] = complex(float(Wcorr[1+2*col_corr]), float(Wcorr[2+2*col_corr]))
        
        for BinSize in range(1, Nconf):
            if (Nconf % BinSize != 0):
                continue
            
            Nbin = Nconf / BinSize
            
            BinCorr_p0 = [0.0] * Nbin
            makeJK(Corr_p0, BinCorr_p0, BinSize)
            
            if (CorEffFlg == 1):
                BinCorr_p1 = [0.0] * Nbin
                makeJK(Corr_p1, BinCorr_p1, BinSize)
                
                BinEffmass = [0.0] * Nbin
                for ibin in range(Nbin):
                    BinEffmass[ibin] = - cmath.log(BinCorr_p1[ibin] / BinCorr_p0[ibin])
            
            if (CorEffFlg == 0):
                makeJKerr(BinCorr_p0, BinSize)
            else:
                makeJKerr(BinEffmass, BinSize)
        print
    quit()
