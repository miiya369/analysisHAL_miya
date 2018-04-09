#!/usr/bin/python

if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 3):
        print "usage: %s [Tshift] [Ens]" % argv[0]
        quit()
    
    Tshift = int(argv[1])
    Ens    = int(argv[2])
    
    OrgParam = "/bghome/scnfqcd/tkymiya/work_Bridge1.3+HF/pacs-cs_ens%d.RHQ/Param/quark_Tshift000.yaml" % Ens
    
    for line in open(OrgParam, 'r'):
        words = line.split()
        
        if (len(words) > 0):
            if (words[0] == "source_position"):
                print "  source_position : [0,0,0,%d]" % Tshift
            else:
                print line.strip("\n")
        else:
            print line.strip("\n")
