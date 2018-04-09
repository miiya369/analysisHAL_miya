#!/usr/bin/python

def elapse_time_sec(hms1, hms2):
    hms1_sp = hms1.split(":")
    hms2_sp = hms2.split(":")
    
    sec1 = int(hms1_sp[0])*60*60 + int(hms1_sp[1])*60 + int(hms1_sp[2])
    sec2 = int(hms2_sp[0])*60*60 + int(hms2_sp[1])*60 + int(hms2_sp[2])
    
    elapse_sec = sec2 - sec1
    
    return elapse_sec

def convert_sec_hms(second):
    diff_sec  = (second % 3600) % 60
    diff_min  = ((second-diff_sec)/60) % 60
    diff_hour = (second-diff_min*60-diff_sec) / 3600
    
    elapse_time_ret = "%d:%d:%d" % (diff_hour, diff_min, diff_sec)
    
    return elapse_time_ret


if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 2):
        print "usage: %s [ifile]" % argv[0]
        quit()
    
    ifile = open(argv[1], 'r')
    lines = ifile.readlines()
    ifile.close()
    
    time  = [0, 0, 0]
    quark = ("ud quark", "s  quark", "c  quark")
    count = [0, 0, 0]
    
    for line in lines:
        words = line.split()
        if (len(words) > 0):
            if (words[0] == "@@@"):
                if (words[1] == "solver(start):"):
                    if (words[8] == "ud"):
                        time1 = words[5]
                        quark_idx = 0
                        count[quark_idx] += 1
                    elif (words[8] == "s"):
                        time1 = words[5]
                        quark_idx = 1
                        count[quark_idx] += 1
                    elif (words[8] == "c"):
                        time1 = words[5]
                        quark_idx = 2
                        count[quark_idx] += 1
                if (words[1] == "solver(end):"):
                    time2 = words[5]
                    time[quark_idx] += elapse_time_sec(time1, time2)
    
    for i in range(0, 3):
        if (count[i] != 0):
            time[i] /= count[i]
            print "%s %s" % (quark[i], convert_sec_hms(time[i]))
