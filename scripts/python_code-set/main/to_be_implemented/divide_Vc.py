#!/usr/bin/python                                                                                                      

from __future__ import print_function

import numpy as np
from struct import pack, unpack

### =========================== Main =========================== ###

### Choose the data which will be averaged
def main (ofpath, ifpaths):
    #return average_textdata(ofpath, ifpaths)
    #return average_bindata (ofpath, ifpaths)
    #return average_wave48  (ofpath, ifpaths)
    return -1

### ============================================================ ###

def average_bindata (ofpath, ifpaths):
    ave_data = np.mean(np.array([np.fromfile(ifpaths[i], '>d') for i in range(len(ifpaths))]), axis=0)
    with open(ofpath, 'wb') as ofdata:
        ofdata.write(pack('>%dd' % len(ave_data), *ave_data))
    return 0

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    from sys     import exit, argv; argc = len(argv)
    from os.path import basename
    
    if (argc != 5):
        exit("usage: python %s [ifname: Vc 1S0] [ifname: Vc 3S1] [ofname: V0] [ofname: Vsigma]" % basename(argv[0]))
    
    ifpath_S0 = argv[1].strip()
    ifpath_S1 = argv[2].strip()
    ofpath_V0 = argv[3].strip()
    ofpath_Vs = argv[4].strip()
    
    if (main(ifpath_S0, ifpath_S1, ofpath_V0, ofpath_Vs) != 0):
        exit("ERROR EXIT.")




if [ $# -ne 4 ]; then
    echo "usage: sh `basename $0` [ifile: Vc 1S0] [ifile: Vc 3S1] [ofile: V0] [ofile: Vsigma]"
    exit 1
fi

function calc() {
    awk "BEGIN {print $1}"
    return $?
}

if_1S=$1
if_3S=$2
of_V0=$3
of_Vs=$4

add_wave=$HOME/python_code-set/main_sub/main_Add.BinData.py

$add_wave $of_V0 `calc +3/4` $if_3S `calc +1/4` $if_1S
$add_wave $of_Vs `calc +1/4` $if_3S `calc -1/4` $if_1S
