#!/usr/bin/python

# author: Takaya Miyamoto (YITP, Kyoto Univ.)
# since : Tue Jul 26 22:58:52 BST 2016
# brief : Make the Jack Knife Samples

######### CHANGE HERE #########
txyzSHIFT = "000.000.000.000"
OpType    = "PP"
SnkRela   = "CG05"
SrcRela   = "CG05"

Endian = 0 # 0 = Little, 1 = Big
###############################

Fsize    = 0
WaveSize = 0
Nconf    = 0

def SetPath(CONF, TIME):
    return "%s/NBSwave.+%03d+%s.%s.%s_%s.%s_%s" \
        % (CONF, TIME, txyzSHIFT, CONF, OpType, SnkRela, OpType, SrcRela),

def SetSize(Fname):
    global Fsize, WaveSize
    
    Fsize = os.path.getsize(Fname)
    if (Fsize % 8 != 0):
        sys.exit("Invalid file size: fsize % sizeof(double) != 0")
    WaveSize = Fsize / 8

def ReadWaveAll(Fname, WAVE, CONF_NUM):
    fsize = os.path.getsize(Fname)    
    if (Fsize != fsize):
        sys.exit("ERROR: Defferent file size, %s" % Fname)
    
    ifwave = open(Fname, 'rb')
    for n in range(WaveSize):
        if (Endian == 0):
            WAVE[n+Nconf*CONF_NUM] = struct.unpack('<d', ifwave.read(8))[0]
        elif (Endian == 1):
            WAVE[n+Nconf*CONF_NUM] = struct.unpack('>d', ifwave.read(8))[0]
        else:
            sys.exit("ERROR: Invalid endian index, %d" % Endian)
    ifwave.close()

def WriteWave(Fname, WAVE):
    ofwave = open(Fname, 'wb')
    for n in range(WaveSize):
        if (Endian == 0):
            ofwave.write(struct.pack('<d', WAVE[n]))
        elif (Endian == 1):
            ofwave.write(struct.pack('>d', WAVE[n]))
        else:
            sys.exit("ERROR: Invalid endian index, %d" % Endian)
    ofwave.close()

if __name__ == "__main__":
    import sys
    import math
    import os.path
    import struct
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 7):
        print("usage: %s [conf list] [bin size] [JKconf base] [idir] [odir] [time slice]"
              % argv[0])
        sys.exit(1)
    
    Flist    = open(argv[1].strip())
    ConfList = Flist.readlines()
    Nconf    = len(ConfList)
    Flist.close()
    
    binSize  = int(argv[2].strip())

    if (Nconf % binSize != 0):
        sys.exit("ERROR: Invalid #.conf and bin size, %d % $d != 0" % (Nconf, binSize))
    Nbin = Nconf / binSize
    
    obase = argv[3].strip()
    idir  = argv[4].strip()
    odir  = argv[5].strip()
    
    TimeSlice = int(argv[6].strip())
    
    SetSize(idir+"/"+SetPath(ConfList[0].strip(), TimeSlice)[0])
    Wave_all = [0.0] * (WaveSize*Nconf)
    Wave_sum = [0.0] *  WaveSize
    Wave_tmp = [0.0] *  WaveSize

    print("#.conf = %d, bin size = %d, #.bin = %d, time = %d"
          % (Nconf, binSize, Nbin, TimeSlice))
    print("input : %s" % idir)
    print("output: %s" % odir)
    
    for i in range(Nconf):
        ReadWaveAll(idir+"/"+SetPath(ConfList[i].strip(), TimeSlice)[0], Wave_all, i)
        for n in range(WaveSize):
            Wave_sum[n] += Wave_all[n+Nconf*i]
    
    for i in range(Nbin):
        for n in range(WaveSize):
            Wave_tmp[n] = Wave_sum[n]
        
        for n in range(WaveSize):
            for j in range(binSize):
                Wave_tmp[n] -= Wave_all[n+Nconf*(j+binSize*i)]
            Wave_tmp[n] /= float(Nconf - binSize)
        
        WriteWave(odir+"/"+SetPath(obase+("_bin%03d-%03d"%(binSize, i)), TimeSlice)[0]
                  , Wave_tmp)
