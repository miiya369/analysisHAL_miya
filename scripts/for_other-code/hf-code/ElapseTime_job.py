#!/usr/bin/python

def calc_time_sec(hms1, ope, hms2):
    hms1_sp = hms1.split(":")
    hms2_sp = hms2.split(":")
    
    sec1 = int(hms1_sp[0])*60*60 + int(hms1_sp[1])*60 + int(hms1_sp[2])
    sec2 = int(hms2_sp[0])*60*60 + int(hms2_sp[1])*60 + int(hms2_sp[2])
    
    if (ope == "+"):
        elapse_sec = sec1 + sec2
    elif (ope == "-"):
        elapse_sec = sec1 - sec2
    
    return elapse_sec

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
        print "usage: %s [ifile1] [ifile2] ..." % argv[0]
        quit()
    
    count = 0
    check = 0
    mean  = 0
    mean2 = 0
    
    for i in range(argc):
        ifile = open(argv[i], 'r')
        lines = ifile.readlines()
        ifile.close()
        check = 0
        
        for line in lines:
            words = line.split()
            if (len(words) > 0):
                if(words[0] == "@@@"):
                    if(words[1] == "JOB(START):"):
                        job_s_day = words[4]
                        job_s     = words[5]
                        check    += 1
                    if(words[1] == "JOB(END),"):
                        job_e_day = words[4]
                        job_e     = words[5]
                        check    += 1
            
        if (check == 2):
            if (job_s_day == job_e_day):
                elapse = calc_time_sec(job_e, "-", job_s)
            else:
                tmp    = calc_time_sec("24:00:00", "-", job_s)
                tmp2   = convert_sec_hms(tmp)
                elapse = calc_time_sec(job_e, "+", tmp2)                
            
            mean  += elapse
            mean2 += elapse * elapse
            count += 1

                    
    if (count != 0):
        mean  /= count
        mean2 /= count
    
    err = 0
    if (count != 1):
        err = int( math.sqrt((mean2 - mean*mean) / (count-1)) )
    
    print "Elapse time of JOB: %s +/- %s (%d)" \
        % (convert_sec_hms(mean), convert_sec_hms(err), count)
