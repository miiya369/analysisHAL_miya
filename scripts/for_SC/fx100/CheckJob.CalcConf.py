#!/usr/bin/python

if __name__ == "__main__":
    import sys
    import math
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc == 1):
        print "usage: %s [job.out file] ..." % argv[0]
        quit()
    
    for i in range(argc-1):
        txyz = ""
        srcT = 0
        tmps = ""
        
        flg1 = 0
        flg2 = 0
        
        ConfNum = 0
        
        ifile = open(argv[i+1], 'r')
        lines = ifile.readlines()
        ifile.close()
        
        for line in lines:
            words = line.split()
            if (len(words) > 0):
                if (words[0] == "Hadron"):
                    if (words[1] == "Save"):
                        if (words[2] == "path"):
                            tmps = words[3].split(".")
                            txyz = tmps[5].split("_")[1]
                            srcT = int(tmps[4].split("_")[1])
            
                if (words[0] == "traj="):
                    if (flg1 == 0):
                        flg1    = 1
                        ConfNum = int(words[1])
                    else:
                        print "%03d %d %s Failed" % (srcT, ConfNum, txyz)
                        ConfNum = int(words[1])
                
                if (words[0] == "NBS_4Boct::~NBS_4Boct():"):
                    flg2 = 1
                
                if (words[0] == "Hadron::run_2pt():"):
                    if (words[2] == "correlator(end):"):
                        if (flg2 == 1):
                            flg2 = 0
                            if (flg1 == 1):
                                print "%03d %d %s Successful" % (srcT, ConfNum, txyz)
                                flg1 = 0
        if (flg1 == 1):
            print "%03d %d %s Failed" % (srcT, ConfNum, txyz)
