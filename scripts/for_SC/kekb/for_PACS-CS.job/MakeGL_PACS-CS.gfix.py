#!/usr/bin/python

if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 4):
        print "usage: %s [Conf Start] [Conf End] [Ens]" % argv[0]
        quit()

    ENS = argv[3]
    
    if (ENS == "1"):
        GFdir = "/bgdata4/scnfqcd/miyamoto/conf/pacs-cs/Kud01370000Ks01364000.gfix"
    elif (ENS == "2"):
        GFdir = "/bgdata4/scnfqcd/miyamoto/conf/pacs-cs/Kud01372700Ks01364000.gfix"  
    elif (ENS == "3"):
        GFdir = "/bgdata4/scnfqcd/miyamoto/conf/pacs-cs/Kud01375400Ks01364000.gfix"
    else:
        quit()
        
    GFlst = "/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/pacs-cs_ens%s.gfix.conf.lst" % ENS
    
    ifile = open(GFlst, 'r')
    lines = ifile.readlines()
    ifile.close()
    
    Nline = len(lines)
    
    for l in range(int(argv[1]), int(argv[2])+1, 1):
        print "%s/%s" % (GFdir, lines[l-1].strip())
