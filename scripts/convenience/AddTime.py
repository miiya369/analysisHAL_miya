#!/usr/bin/python

def elapse_time_sec(hms1, hms2):
    hms1_sp = hms1.split(":")
    hms2_sp = hms2.split(":")
    
    sec1 = int(hms1_sp[0])*60*60 + int(hms1_sp[1])*60 + int(hms1_sp[2])
    sec2 = int(hms2_sp[0])*60*60 + int(hms2_sp[1])*60 + int(hms2_sp[2])
    
    elapse_sec = sec2 - sec1
    
    return elapse_sec

def convert_hms_sec(hms):
    hms_sp = hms.split(":")
    
    sec = int(hms_sp[0])*60*60 + int(hms_sp[1])*60 + int(hms_sp[2])
    
    return sec

def convert_sec_hms(second):
    diff_sec  = (second % 3600) % 60
    diff_min  = ((second-diff_sec)/60) % 60
    diff_hour = (second-diff_min*60-diff_sec) / 3600
    
    elapse_time_ret = "%d:%d:%d" % (diff_hour, diff_min, diff_sec)
    
    return elapse_time_ret


if __name__ == "__main__":
    import sys
    import math
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc == 1):
        print "usage: %s [hms1] [hms2] ..." % argv[0]
        quit()
    
    time = 0
    
    for i in range(argc-1):
        time += convert_hms_sec(argv[i+1])

    print "%s" % convert_sec_hms(time)
