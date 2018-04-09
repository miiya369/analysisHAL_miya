#!/usr/bin/python

### Auther: Takaya Miyamoto
### Date  : Mon Oct 10 18:19:30 JST 2016
###
### Brief : Plot the R-correlator 
###       : or effective mass of R-corr
###       :

Nconf = 1
Tsize = 1

col_corr = 0 # <- The column of correlator-file to read
             # (For complex, 0 => column 1 & 2, 1 => column 3 & 4, ... are read)

def makeJK(DATAorg, DATAjk, bin_size):
    if (Nconf != 1 and bin_size != 0):
        Nbin = Nconf / bin_size
        
        for it in range(Tsize):
            TmpSum = complex(0.0, 0.0)
            
            for i in range(Nconf):
                TmpSum += DATAorg[it + Tsize * i]
            
            for i in range(Nbin):
                TmpSum2 = TmpSum
                
                for j in range(bin_size):
                    TmpSum2 -= DATAorg[it + Tsize * (j + bin_size * i)]
                
                DATAjk[it + Tsize * i] = TmpSum2 / float(Nconf-bin_size)
    else:
        for i in range(Tsize * Nconf):
            DATAjk[i] = DATAorg[i]

def makeJKerr(DATA, bin_size, COR_Eff_FLG):
    if (bin_size != 0):
        Nbin = Nconf / bin_size
    else:
        Nbin = Nconf
    
    for it in range(Tsize-COR_Eff_FLG):
        MeanR  = 0.0
        Mean2R = 0.0
        MeanI  = 0.0
        Mean2I = 0.0
        
        for ibin in range(Nbin):
            if (COR_Eff_FLG == 0):
                DATA_tmp = DATA[it + Tsize * ibin]
            else:
                if (DATA[it+1 + Tsize * ibin] == complex(0.0, 0.0) or
                    DATA[it   + Tsize * ibin] == complex(0.0, 0.0)):
                    DATA_tmp = complex(0.0, 0.0)
                else:
                    DATA_tmp = - cmath.log(DATA[it+1 + Tsize * ibin] / DATA[it + Tsize * ibin])
            
            MeanR  += DATA_tmp.real
            Mean2R += DATA_tmp.real ** 2
            MeanI  += DATA_tmp.imag
            Mean2I += DATA_tmp.imag ** 2
        
        MeanR  /= float(Nbin)
        Mean2R /= float(Nbin)
        MeanI  /= float(Nbin)
        Mean2I /= float(Nbin)
        ErrR    = math.sqrt((Nbin-1) * (Mean2R - MeanR ** 2))
        ErrI    = math.sqrt((Nbin-1) * (Mean2I - MeanI ** 2))
        
        print "%d %e %e %e %e" % (it, MeanR, ErrR, MeanI, ErrI)

if (__name__ == "__main__"):
    import sys
    import math
    import cmath
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 7):
        print("usage: %s [conf-list] [corr-file 1] [corr-file 2] [corr-multi-file] [BinSize (0=Use org-data)] [0=OrgCorr, 1=Effmass]" % argv[0])
        print("     : At the corr-file, \"REPCONF\" will be replaced by conf-name in the conf-list.")
        quit()
    
    ConfList  = argv[1]
    Corr1File = argv[2]
    Corr2File = argv[3]
    CorrMFile = argv[4]
    BinSize   = int(argv[5])
    CorEffFlg = int(argv[6])
    
    if (CorEffFlg != 0 and CorEffFlg != 1):
        print("usage: %s [conf-list] [corr-file 1] [corr-file 2] [corr-multi-file] [BinSize (0=Use org-data)] [0=OrgCorr, 1=Effmass]" % argv[0])
        print("     : At the corr-file, \"REPCONF\" will be replaced by conf-name in the conf-list.")
        quit()
    
    fconf = open(ConfList.strip(), 'r')
    Confs = fconf.readlines()
    Nconf = len(Confs)
    fconf.close()
    
    if (BinSize != 0 and Nconf % BinSize != 0):
        print("#.conf is indivisible by bin size.")
        quit()
    
    fcorrM = open(CorrMFile.replace('REPCONF', Confs[0].strip()), 'r')
    LcorrM = fcorrM.readlines()
    Tsize = len(LcorrM)
    fcorrM.close()
    
    print("# Setup:")
    print("#    - Conf list   = %s" % ConfList)
    print("#    - Corr 1 file = %s" % Corr1File)
    print("#    - Corr 2 file = %s" % Corr2File)
    print("#    - Corr M file = %s" % CorrMFile)
    print("#    - #.conf      = %d" % Nconf)
    print("#    - time size   = %d" % Tsize)
    if (BinSize != 0):
        print("#    - BinSize   = %d" % BinSize)
    else:
        print("#    - BinSize   = Use org-data")
    if (CorEffFlg == 0):
        print("#    - Plot for  : correlator")
    else:
        print("#    - Plot for  : effective mass")
    print("#")
    
    Corr1 = [0.0] * Tsize * Nconf
    Corr2 = [0.0] * Tsize * Nconf
    CorrM = [0.0] * Tsize * Nconf
    
    for iconf in range(Nconf):
        fcorr1 = open(Corr1File.replace('REPCONF', Confs[iconf].strip()), 'r')
        fcorr2 = open(Corr2File.replace('REPCONF', Confs[iconf].strip()), 'r')
        fcorrM = open(CorrMFile.replace('REPCONF', Confs[iconf].strip()), 'r')
        
        Lcorr1 = fcorr1.readlines()
        Lcorr2 = fcorr2.readlines()
        LcorrM = fcorrM.readlines()
        
        fcorr1.close()
        fcorr2.close()
        fcorrM.close()
        
        for it in range(Tsize):
            Wcorr1 = Lcorr1[it].split()
            Wcorr2 = Lcorr2[it].split()
            WcorrM = LcorrM[it].split()
            
            Corr1[it + Tsize * iconf] = complex(float(Wcorr1[1+2*col_corr]), float(Wcorr1[2+2*col_corr]))
            Corr2[it + Tsize * iconf] = complex(float(Wcorr2[1+2*col_corr]), float(Wcorr2[2+2*col_corr]))
            CorrM[it + Tsize * iconf] = complex(float(WcorrM[1+2*col_corr]), float(WcorrM[2+2*col_corr]))
    
    if (BinSize != 0):
        Nbin = Nconf / BinSize
    else:
        Nbin = Nconf
    
    JKcorr1 = [0.0] * Tsize * Nbin
    JKcorr2 = [0.0] * Tsize * Nbin
    JKcorrM = [0.0] * Tsize * Nbin
    
    makeJK(Corr1, JKcorr1, BinSize)
    makeJK(Corr2, JKcorr2, BinSize)
    makeJK(CorrM, JKcorrM, BinSize)
    
    Rcorr = [0.0] * Tsize * Nbin
    
    for ibin in range(Nbin):
        for it in range(Tsize):
            Rcorr[it + Tsize * ibin] = (JKcorrM[it + Tsize * ibin] / 
                                        (JKcorr1[it + Tsize * ibin] * JKcorr2[it + Tsize * ibin]))
    
    makeJKerr(Rcorr, BinSize, CorEffFlg)
