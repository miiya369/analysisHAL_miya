#!/usr/bin/python
# -*- coding: utf-8 -*-

### Auther: Takaya Miyamoto
### Date  : Sun May 28 00:27:38 JST 2017
### Brief : Plot the dispersion relation

### !!! WARNING !!!
### Forward / Backward average for Baryon corr.
### has not implemented yet.

import sys, os
import numpy as np
import subprocess as sp
import math

from scipy.optimize import curve_fit

iDir     = None

had_name = 'pion'
t_min    = 10
t_max    = 20
Lsize    = 32
Do_JK    = False
Ave_FB   = False
bin_size = 1

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [path to correlator.dir] {option}\n" % ARGV[0])
    print("option:")
    print("      -had      [The name of hadron to calculate  ] Default ="),; print had_name
    print("      -t_min    [The minimum range of time for fit] Default ="),; print t_min
    print("      -t_max    [The maximum range of time for fit] Default ="),; print t_max
    print("      -lsize    [The lattice size (Space)         ] Default ="),; print Lsize
    print("      -do_JK     Make the Jack-Knife sample in the code."     ),; print
    print("      -ave_fb    Read the anti-baryon and take average."      ),; print
    print("      -bin_size [Bin Size]                          Default ="),; print bin_size
    print; quit()

def set_args(ARGC, ARGV):
    global iDir, t_min, t_max, Lsize, had_name, Do_JK, Ave_FB, bin_size
    
    if (ARGV[1][0] == '-'):
        usage(ARGV)
    
    iDir = ARGV[1].strip()
    
    for i in range(2, ARGC):
        if (ARGV[i][0] == '-'):
            if   (ARGV[i] == '-had'):
                had_name = ARGV[i+1].strip()
            elif (ARGV[i] == '-t_min'):
                t_min = int(ARGV[i+1])
            elif (ARGV[i] == '-t_max'):
                t_max = int(ARGV[i+1])
            elif (ARGV[i] == '-lsize'):
                Lsize = int(ARGV[i+1])
            elif (ARGV[i] == '-do_JK'):
                Do_JK = True
            elif (ARGV[i] == '-ave_fb'):
                Ave_FB = True
            elif (ARGV[i] == '-bin_size'):
                bin_size = int(ARGV[i+1])
            else:
                print("\nERROR: Invalid argument '%s'" % ARGV[i])
                usage(ARGV)

def check_args():
    print("# === Check Arguments ===")
    print("# idir      ="),; print iDir
    print("# fit t_min ="),; print t_min
    print("# fit t_max ="),; print t_max
    print("# lat.size  ="),; print Lsize
    print("# had name  ="),; print had_name
    print("# Do JK     ="),; print Do_JK
    print("# Ave. f/b  ="),; print Ave_FB
    print("# bin size  ="),; print bin_size
    print("# =======================")

###### Main part
if (__name__ == "__main__"):
    argv = sys.argv; argc = len(argv)
    
    if (argc < 2):
        usage(argv)
    
    set_args(argc, argv); check_args(); #quit()
    
### Input Correlator files ###
    DataList = sp.check_output(['find', iDir, '-type', 'f', '-name', had_name+'*']).split()
    Ndata    = len(DataList)
    
    if (Ave_FB):
        DataList_bw = sp.check_output(['find', iDir, '-type', 'f', '-name', "anti"+had_name+'*']).split()
        if (len(DataList_bw) != Ndata):
            print("\nERROR: #.Data_fw and #.Data_bw are different, exit.\n"); quit()
    
    tmpFile = open(DataList[0].strip(), 'r')
    tmpData = tmpFile.readlines(); tmpFile.close()
    Tsize   = len(tmpData); Ncolm = len(tmpData[0].split())
    del tmpData
    
    print("# N.data    ="),; print Ndata
    print("# T_size    ="),; print Tsize
    print("# N.column  ="),; print Ncolm
    
    if (Ncolm != 11):
        print("\nERROR: Invalid #.colomn to read: #.colomn must be 11, exit.\n"); quit()
        
    Corr = np.empty((Tsize, 5, Ndata))
    
    for idata in range(Ndata):
        tmpFile = open(DataList[idata].strip(), 'r')
        tmpData = tmpFile.readlines(); tmpFile.close()
        print("# Read: %s" % DataList[idata].strip())
        if (len(tmpData) != Tsize):
            print("\nERROR: T_size is differ in '%s', exit.\n" % DataList[idata].strip()); quit()
        
        if (Ave_FB):
            tmpFile    = open(DataList_bw[idata].strip(), 'r')
            tmpData_bw = tmpFile.readlines(); tmpFile.close()
            print("# Read: %s" % DataList_bw[idata].strip())
            if (len(tmpData_bw) != Tsize):
                print("\nERROR: T_size is differ in '%s', exit.\n" % DataList_bw[idata].strip()); quit()
            
            for it in range(Tsize):
                tmpColm_fw = tmpData   [       it         ].split()
                tmpColm_bw = tmpData_bw[(Tsize-it) % Tsize].split()
                if (len(tmpColm_fw) != 11 or len(tmpColm_bw) != 11):
                    print("\nERROR: Invalid #.colomn to read: #.colomn must be 11, exit.\n"); quit()
                for imom in range(5):
                    Corr[it, imom, idata] = (float(tmpColm_fw[imom * 2 + 1]) +
                                             float(tmpColm_bw[imom * 2 + 1])) / 2.0
            del tmpData
            del tmpData_bw
        else:
            for it in range(Tsize):
                tmpColm = tmpData[it].split()
                if (len(tmpColm) != 11):
                    print("\nERROR: Invalid #.colomn to read: #.colomn must be 11, exit.\n"); quit()
                for imom in range(5):
                    Corr[it, imom, idata] = float(tmpColm[imom * 2 + 1])
            del tmpData
    #quit()
    
### Construct the Jack-Knife sample ###
    if (Do_JK):
        if (Ndata % bin_size != 0):
            print("\nERROR: #.data cannot be divided by bin-size: %d %% %d != 0, exit.\n" % 
                  (Ndata, bin_size)); quit()
        
        Ndata_org = Ndata; Ndata /= bin_size
        Corr_org  = Corr; del Corr
        Corr      = np.empty((Tsize, 5, Ndata))
        
        for it in range(Tsize):
            for imom in range(5):
                tmp1_d = np.sum(Corr_org[it, imom, :])
                
                for ibin in range(Ndata):
                    tmp2_d = tmp1_d
                    
                    for ibin_size in range(bin_size):
                        tmp2_d -= Corr_org[it, imom, (ibin_size + bin_size * ibin)]
                        
                    Corr[it, imom, ibin] = tmp2_d / float(Ndata_org-bin_size)
        
        print("# Jack-Knife sample are generated: #.data = %d --> %d" % (Ndata_org, Ndata))
    #quit()
    
### Calculate the effective mass by fitting ###
    Effmass = np.empty((5, Ndata))
    tmpErr  = np.empty((t_max-t_min+1))
    
    for imom in range(5):
        for it in range(t_min, t_max+1):
            tmpErr[it-t_min] = np.sqrt(abs(float(Ndata-1) * 
                                           (np.mean(Corr[it, imom, :]**2) - 
                                            np.mean(Corr[it, imom, :]   )**2)
                                           ))
        for idata in range(Ndata):
            tmpParams, cov = curve_fit(lambda t,a,b: a*np.exp(-b*t), 
                                       np.array(range(t_min, t_max+1)), 
                                       Corr[t_min:t_max+1, imom, idata], sigma=tmpErr)
            
            Effmass[imom, idata] = tmpParams[1]
    
### Calculate the dispersion relations and output ###
    disp_rela = np.empty(Ndata)
    lat_mom   = np.array([0, 1.0, np.sqrt(2.0), np.sqrt(3.0), 2.0]) * 2.0*np.pi / Lsize
    print("#\n# p^2 [Latt.Unit], E^2(p)-E^2(0) [Latt.Unit] {mean, error}")
    for imom in range(1, 5):
        for idata in range(Ndata):
            disp_rela[idata] = Effmass[imom, idata]**2 - Effmass[0, idata]**2
        
        mean = np.mean(disp_rela)
        err  = np.sqrt(abs(float(Ndata-1) * (np.mean(disp_rela**2)-np.mean(disp_rela)**2)))
        print("%f %1.16e %1.16e" % (lat_mom[imom]**2, mean, err))
