#!/usr/bin/python

# Author: Takaya Miyamoto
# Date  : Fri Sep 16 16:28:17 JST 2016
# Brief : Compare the column of two files and show difference

if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 4):
        print "usage: %s [Column number] [file1] [file2]" % argv[0]
        sys.exit(1)
    
    Col = int(argv[1])
    
    if (Col < 1):
        print "Column number shold be set greater than 0"
        sys.exit(1)
    
    file1  = open(argv[2], 'r')
    lines1 = file1.readlines()
    file1.close()
    
    file2  = open(argv[3], 'r')
    lines2 = file2.readlines()
    file2.close()
    
    len1 = len(lines1)
    len2 = len(lines2)
    
    if (len1 != len2):
        print "#.line of each file is different."
        sys.exit(1)
    
    for i in range(len1):
        str1 = lines1[i].split()
        str2 = lines2[i].split()
        
        if (str1[Col-1] != str2[Col-1]):
            print "%s and %s are differ" % (argv[2], argv[3])
            sys.exit(1)

    sys.exit(0)
