#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Tue May 30 19:15:14 JST 2017
### Brief : Calculate the NLO potentials for NxN coupled channel.

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np

from MiscFuncs.Misc                import *
from MiscFuncs.DataStatistics      import make_mean_err
from MiscFuncs.IO_BinData          import input_bin_data, output_bin_data
from PotentialNLO.CalcPotentialNLO import calc_potential_NLO

oFbase = None

###### Main part
if __name__ == "__main__":
    argv = sys.argv; argc = len(argv)
    if (argc != 4):
        print("usage: %s [idir] [time] [odir]" % argv[0]); quit()
    
    fbase  = argv[1].strip()
    time   = int(argv[2].strip())
    oFbase = argv[3].strip()
    
    # ifiles[time, pot/Rcorr/LapR, exp/wall, channel(0-3)]
    
    ifiles=np.array(
        [
            [
                [
                    (fbase+"/LN-LN_pot_expo_t%d"%(time-1), fbase+"/LN-SN_pot_expo_t%d"%(time-1), 
                     fbase+"/SN-LN_pot_expo_t%d"%(time-1), fbase+"/SN-SN_pot_expo_t%d"%(time-1)),
                    (fbase+"/LN-LN_pot_wall_t%d"%(time-1), fbase+"/LN-SN_pot_wall_t%d"%(time-1),
                     fbase+"/SN-LN_pot_wall_t%d"%(time-1), fbase+"/SN-SN_pot_wall_t%d"%(time-1))
                    ], [
                    (fbase+"/LN-LN_Rcr_expo_t%d"%(time-1), fbase+"/LN-SN_Rcr_expo_t%d"%(time-1),
                     fbase+"/SN-LN_Rcr_expo_t%d"%(time-1), fbase+"/SN-SN_Rcr_expo_t%d"%(time-1)),
                    (fbase+"/LN-LN_Rcr_wall_t%d"%(time-1), fbase+"/LN-SN_Rcr_wall_t%d"%(time-1),
                     fbase+"/SN-LN_Rcr_wall_t%d"%(time-1), fbase+"/SN-SN_Rcr_wall_t%d"%(time-1))
                    ], [
                    (fbase+"/LN-LN_LaR_expo_t%d"%(time-1), fbase+"/LN-SN_LaR_expo_t%d"%(time-1),
                     fbase+"/SN-LN_LaR_expo_t%d"%(time-1), fbase+"/SN-SN_LaR_expo_t%d"%(time-1)),
                    (fbase+"/LN-LN_LaR_wall_t%d"%(time-1), fbase+"/LN-SN_LaR_wall_t%d"%(time-1),
                     fbase+"/SN-LN_LaR_wall_t%d"%(time-1), fbase+"/SN-SN_LaR_wall_t%d"%(time-1))
                    ]
                ],
            [
                [
                    (fbase+"/LN-LN_pot_expo_t%d"%(time), fbase+"/LN-SN_pot_expo_t%d"%(time),
                     fbase+"/SN-LN_pot_expo_t%d"%(time), fbase+"/SN-SN_pot_expo_t%d"%(time)),
                    (fbase+"/LN-LN_pot_wall_t%d"%(time), fbase+"/LN-SN_pot_wall_t%d"%(time),
                     fbase+"/SN-LN_pot_wall_t%d"%(time), fbase+"/SN-SN_pot_wall_t%d"%(time))
                    ], [
                    (fbase+"/LN-LN_Rcr_expo_t%d"%(time), fbase+"/LN-SN_Rcr_expo_t%d"%(time),
                     fbase+"/SN-LN_Rcr_expo_t%d"%(time), fbase+"/SN-SN_Rcr_expo_t%d"%(time)),
                    (fbase+"/LN-LN_Rcr_wall_t%d"%(time), fbase+"/LN-SN_Rcr_wall_t%d"%(time),
                     fbase+"/SN-LN_Rcr_wall_t%d"%(time), fbase+"/SN-SN_Rcr_wall_t%d"%(time))
                    ], [
                    (fbase+"/LN-LN_LaR_expo_t%d"%(time), fbase+"/LN-SN_LaR_expo_t%d"%(time),
                     fbase+"/SN-LN_LaR_expo_t%d"%(time), fbase+"/SN-SN_LaR_expo_t%d"%(time)),
                    (fbase+"/LN-LN_LaR_wall_t%d"%(time), fbase+"/LN-SN_LaR_wall_t%d"%(time),
                     fbase+"/SN-LN_LaR_wall_t%d"%(time), fbase+"/SN-SN_LaR_wall_t%d"%(time))
                    ]
                ],
            [
                [
                    (fbase+"/LN-LN_pot_expo_t%d"%(time+1), fbase+"/LN-SN_pot_expo_t%d"%(time+1),
                     fbase+"/SN-LN_pot_expo_t%d"%(time+1), fbase+"/SN-SN_pot_expo_t%d"%(time+1)),
                    (fbase+"/LN-LN_pot_wall_t%d"%(time+1), fbase+"/LN-SN_pot_wall_t%d"%(time+1),
                     fbase+"/SN-LN_pot_wall_t%d"%(time+1), fbase+"/SN-SN_pot_wall_t%d"%(time+1))
                    ], [
                    (fbase+"/LN-LN_Rcr_expo_t%d"%(time+1), fbase+"/LN-SN_Rcr_expo_t%d"%(time+1),
                     fbase+"/SN-LN_Rcr_expo_t%d"%(time+1), fbase+"/SN-SN_Rcr_expo_t%d"%(time+1)),
                    (fbase+"/LN-LN_Rcr_wall_t%d"%(time+1), fbase+"/LN-SN_Rcr_wall_t%d"%(time+1),
                     fbase+"/SN-LN_Rcr_wall_t%d"%(time+1), fbase+"/SN-SN_Rcr_wall_t%d"%(time+1))
                    ], [
                    (fbase+"/LN-LN_LaR_expo_t%d"%(time+1), fbase+"/LN-SN_LaR_expo_t%d"%(time+1),
                     fbase+"/SN-LN_LaR_expo_t%d"%(time+1), fbase+"/SN-SN_LaR_expo_t%d"%(time+1)),
                    (fbase+"/LN-LN_LaR_wall_t%d"%(time+1), fbase+"/LN-SN_LaR_wall_t%d"%(time+1),
                     fbase+"/SN-LN_LaR_wall_t%d"%(time+1), fbase+"/SN-SN_LaR_wall_t%d"%(time+1))
                    ]
                ]
            ]
        )
    
    Nch  = len(ifiles[0,0,0,:]) / 2
    Nsrc = len(ifiles[0,0,:,0])
    Nt   = len(ifiles[:,0,0,0])
    
    Neq  = Nt * Nsrc
    
    tmp_yData, xData, dummy = input_bin_data(ifiles[0, 0, 0, 0])
    
    Ndata = len(tmp_yData[:, 0])
    Nconf = len(tmp_yData[0, :])
    
    del tmp_yData
    del dummy
    
    print("#")
    print("# N.ch   = "),; print Nch
    print("# N.src  = "),; print Nsrc
    print("# N.t    = "),; print Nt
    print("# N.eq   = "),; print Neq
    print("# N.data = "),; print Ndata
    print("# N.conf = "),; print Nconf
    print("#")
    
    Pot = np.empty((Nch, Nch, Ndata, Neq, Nconf))
    Rcr = np.empty((Nch, Nch, Ndata, Neq, Nconf))
    LpR = np.empty((Nch, Nch, Ndata, Neq, Nconf))
    
    for it in range(Nt):
        for isrc in range(Nsrc):
            for i in range(Nch):
                for j in range(Nch):
                    tmpPot, xData_Pot, dummy = input_bin_data(ifiles[it, 0, isrc, j + Nch*i])
                    tmpRcr, xData_Rcr, dummy = input_bin_data(ifiles[it, 1, isrc, j + Nch*i])
                    tmpLpR, xData_LpR, dummy = input_bin_data(ifiles[it, 2, isrc, j + Nch*i])
                    
                    if (len(xData_Pot) != Ndata or
                        len(xData_Rcr) != Ndata or
                        len(xData_LpR) != Ndata):
                        print("ERROR: len(xData) is different, exit."); quit()
                    
                    for idata in range(Ndata):
                        if (xData[idata] != xData_Pot[idata] or
                            xData[idata] != xData_Rcr[idata] or
                            xData[idata] != xData_LpR[idata]):
                            print("ERROR: x-data is different, exit."); quit()
                    
                        for iconf in range(Nconf):
                            Pot[i, j, idata, isrc + Nsrc*it, iconf] = tmpPot[idata, iconf]
                            Rcr[i, j, idata, isrc + Nsrc*it, iconf] = tmpRcr[idata, iconf]
                            LpR[i, j, idata, isrc + Nsrc*it, iconf] = tmpLpR[idata, iconf]
    
    Pot__LO = np.empty((Nch, Nch, Ndata, Nconf))
    Pot_NLO = np.empty((Nch, Nch, Ndata, Nconf))
    
    for iconf in range(Nconf):
        tmpPot__LO, tmpPot_NLO = calc_potential_NLO(Pot[:,:,:,:,iconf], 
                                                    Rcr[:,:,:,:,iconf], 
                                                    LpR[:,:,:,:,iconf])
        if (tmpPot__LO is tmpPot_NLO is None):
            quit()
        for idata in range(Ndata):
            for i in range(Nch):
                for j in range(Nch):
                    Pot__LO[i, j, idata, iconf] = tmpPot__LO[i, j, idata]
                    Pot_NLO[i, j, idata, iconf] = tmpPot_NLO[i, j, idata]
    
    output_bin_data(oFbase+"/LN-LN_CCP__LO_t%d"%(time), Pot__LO[0, 0, :, :], xData)
    output_bin_data(oFbase+"/LN-LN_CCP_NLO_t%d"%(time), Pot_NLO[0, 0, :, :], xData)
    output_bin_data(oFbase+"/LN-SN_CCP__LO_t%d"%(time), Pot__LO[0, 1, :, :], xData)
    output_bin_data(oFbase+"/LN-SN_CCP_NLO_t%d"%(time), Pot_NLO[0, 1, :, :], xData)
    output_bin_data(oFbase+"/SN-LN_CCP__LO_t%d"%(time), Pot__LO[1, 0, :, :], xData)
    output_bin_data(oFbase+"/SN-LN_CCP_NLO_t%d"%(time), Pot_NLO[1, 0, :, :], xData)
    output_bin_data(oFbase+"/SN-SN_CCP__LO_t%d"%(time), Pot__LO[1, 1, :, :], xData)
    output_bin_data(oFbase+"/SN-SN_CCP_NLO_t%d"%(time), Pot_NLO[1, 1, :, :], xData)
    
    for i in range(Nch):
        for j in range(Nch):
            for idata in range(Ndata):
                mean0, err0 = make_mean_err(Pot__LO[i, j, idata, :])
                mean1, err1 = make_mean_err(Pot_NLO[i, j, idata, :])
                print("%lf %e %e %e %e" % (xData[idata], mean0, err0, mean1, err1))
            print("\n")
