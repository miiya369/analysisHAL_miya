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
    import math
    
    argv = sys.argv
    argc = len(argv)
    
    if (argc != 2):
        print "usage: %s [ifile]" % argv[0]
        quit()
    
    ifile = open(argv[1], 'r')
    lines = ifile.readlines()
    ifile.close()

    ttime = [0, 0, 0]
    time  = [0, 0, 0]
    quark = ("Solver 1", "Solver 2", "HF  code")
    count = [0, 0, 0]
    tmp   = 0
    tmp_s = -1

    alltime = 0
    jobtime = 0
    
    for line in lines:
        words = line.split()
        if (len(words) > 1):
            if (words[0] == "@@@"):
                if (words[1] == "solver(start):"):
                    time1  = words[5]
                    tmp_s += 1
                if (words[1] == "solver(end):"):
                    time2         = words[5]
                    tmp_t         = elapse_time_sec(time1, time2)
                    count[tmp_s] += 1
                    time [tmp_s] += tmp_t
                    ttime[tmp_s] += tmp_t * tmp_t
                    if (tmp_s == 1):
                        tmp_s = -1
                
                if (words[1] == "JOB(START):"):
                    job_s = words[5]
                if (words[1] == "JOB(END),"):
                    job_e   = words[5]
                    jobtime = elapse_time_sec(job_s, job_e)
            
            if (words[0] == "Hadron::run_2pt():"):
                if (words[2] == "correlator(begin):"):
                    if (tmp == 0):
                        time1 = words[6]
                        tmp   = 1
                    elif (tmp == 1):
                        tmp = 0
                if (words[2] == "correlator(end):"):
                    if (tmp == 0):
                        time2     = words[6]
                        tmp_t     = elapse_time_sec(time1, time2)
                        count[2] += 1
                        time [2] += tmp_t
                        ttime[2] += tmp_t * tmp_t
    
    for i in range(0, 3):
        if (count[i] != 0):
            time[i]  /= count[i]
            ttime[i] /= count[i]
            
            err = 0
            if (count[i] != 1):
                err = int( math.sqrt((ttime[i]-time[i]*time[i]) / (count[i]-1)) )
            
            print "%s %s +/- %s (%d)" \
            % (quark[i], convert_sec_hms(time[i]), convert_sec_hms(err), count[i])
            
            alltime += time[i] * count[i]
            

    print "\nAllTime = %s" % convert_sec_hms(alltime)
    print "JobTime = %s" % convert_sec_hms(jobtime)
