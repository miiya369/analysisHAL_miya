#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from Fitting.IO_Params import input_params

iFname = None
oFname = None

###### Functions for arguments
def usage(ARGV):
    print("\nusage: %s [ifname (text-data)] [ofname (bin-data)]\n" % ARGV[0])
    quit()

def set_args(ARGC, ARGV):
    global iFname, oFname
    
    iFname = ARGV[1].strip()
    oFname = ARGV[2].strip()

###### Main part
if __name__ == "__main__":
    from struct import pack
    
    argv = sys.argv; argc = len(argv)
    
    if (argc != 3):
        usage(argv)
    
    set_args(argc, argv)
        
    FitFunc_name, Params = input_params(iFname)
    if (FitFunc_name is Params is None):
        quit()
    
    if (FitFunc_name == '1E'):
        FitFunc_num = 2
    elif (FitFunc_name == '1G'):
        FitFunc_num = 3
    elif (FitFunc_name == '2G'):
        FitFunc_num = 4
    elif (FitFunc_name == '3G'):
        FitFunc_num = 5
    elif (FitFunc_name == '1SG'):
        FitFunc_num = 6
    elif (FitFunc_name == '2SG'):
        FitFunc_num = 7
    elif (FitFunc_name == '4G'):
        FitFunc_num = 8
    else:
        print("ERROR: Invalid Fit function name '%s', exit." % FitFunc_name)
        quit()
    
    ofile = open(oFname, 'wb')
    ofile.write(pack('<i', len(Params[0, :])))
    ofile.write(pack('<i', FitFunc_num))
    for iparam in range(len(Params[:, 0])):
        for iconf in range(len(Params[0, :])):
            ofile.write(pack('<d', Params[iparam, iconf]))
    ofile.close()
