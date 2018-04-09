#!/usr/bin/python

import math
import numpy as np

Ncount = 14

enum_totjob = 0
enum_reconf = 1
enum_sol_ud = 2
enum_sol__s = 3
enum_sol__c = 4
enum_2ptcor = 5
enum_2pt_mm = 6
enum_2pt_oo = 7
enum_2pt_dd = 8
enum_2pt_od = 9
enum_nbs_mm = 10
enum_nbs_oo = 11
enum_nbs_dd = 12
enum_nbs_od = 13

enum_name = ("Total JOB",
             "Read Conf",
             "Solver ud",
             "Solver  s",
             "Solver  c",
             "2pt  corr",
             "2pt Me-Me",
             "2pt Oc-Oc",
             "2pt De-De",
             "2pt Oc-De",
             "NBS Me-Me",
             "NBS Oc-Oc",
             "NBS De-De",
             "NBS Oc-De")

# For bridge++ solver
tmp_quark_id = 0

tot_time_sec  = np.zeros(Ncount, dtype=int)
tot_time_sec2 = np.zeros(Ncount, dtype=int)
tmp_time_sec  = np.zeros(Ncount, dtype=int)
num_count     = np.zeros(Ncount, dtype=int)

def main(argc, argv):
    for i in range(argc):
        ifile = open(argv[i].strip(), 'r')
        LINES = ifile.readlines(); ifile.close()
        file_iter(LINES)

def file_iter(lines):
    for line in lines:
        WORDS = line.split()
        if (len(WORDS) > 1):
            if (WORDS[0] == "@@@"):
                count_bridge(WORDS)
            if (WORDS[1] == "@@@"):
                count_hf(WORDS)

def count_bridge(words):
    global tmp_quark_id
    
    if (words[1] == "JOB(START):"):
        tmp_time_sec[enum_totjob] = conv_hmc_to_sec(words[5].strip())
    if (words[1] == "JOB(END),"):
        num_count[enum_totjob] += 1
        add_count(enum_totjob, conv_hmc_to_sec(words[5].strip()))
    
    if (words[1] == "read" and words[2] == "conf(start)"):
        tmp_time_sec[enum_reconf] = conv_hmc_to_sec(words[6].strip())
    if (words[1] == "read" and words[2] == "conf(end)"):
        num_count[enum_reconf] += 1
        add_count(enum_reconf, conv_hmc_to_sec(words[6].strip()))
    
    if (words[1] == "solver(start):"):
        if (words[8] == "ud"):
            tmp_quark_id = 0
        if (words[8] == "s"):
            tmp_quark_id = 1
        if (words[8] == "c"):
            tmp_quark_id = 2
        tmp_time_sec[enum_sol_ud+tmp_quark_id] = conv_hmc_to_sec(words[5].strip())
    if (words[1] == "solver(end):"):
        num_count[enum_sol_ud+tmp_quark_id] += 1
        add_count(enum_sol_ud+tmp_quark_id, conv_hmc_to_sec(words[5].strip()))

def count_hf(words):
    if (words[0] == "Hadron::run_2pt():" and words[2] == "correlator(begin):"):
        tmp_time_sec[enum_2ptcor] = conv_hmc_to_sec(words[6].strip())
    if (words[0] == "Hadron::run_2pt():" and words[2] == "correlator(end):"):
        num_count[enum_2ptcor] += 1
        add_count(enum_2ptcor, conv_hmc_to_sec(words[6].strip()))
    
    if (words[0] == "NBS_2Boct::run_2pt():" and words[2] == "WAVE_2PT[2Boct](start):"):
        tmp_time_sec[enum_2pt_oo] = conv_hmc_to_sec(words[9].strip())
    if (words[0] == "NBS_2Boct::run_2pt():" and words[2] == "WAVE_2PT[2Boct](end):"):
        num_count[enum_2pt_oo] += 1
        add_count(enum_2pt_oo, conv_hmc_to_sec(words[9].strip()))
    
    if (words[0] == "NBS_2Bdow::run_2pt():" and words[2] == "WAVE_2PT[2Bdow](start):"):
        tmp_time_sec[enum_2pt_od] = conv_hmc_to_sec(words[9].strip())
    if (words[0] == "NBS_2Bdow::run_2pt():" and words[2] == "WAVE_2PT[2Bdow](end):"):
        num_count[enum_2pt_od] += 1
        add_count(enum_2pt_od, conv_hmc_to_sec(words[9].strip()))
    
    if (words[0] == "NBS_2Mpsw::run_NBS():" and words[2] == "WAVE[PP](start):"):
        tmp_time_sec[enum_nbs_mm] = conv_hmc_to_sec(words[9].strip())
    if (words[0] == "NBS_2Mpsw::run_NBS():" and words[2] == "WAVE[PP](end):"):
        num_count[enum_nbs_mm] += 1
        add_count(enum_nbs_mm, conv_hmc_to_sec(words[9].strip()))
    
    if (words[0] == "NBS_2Boct::run_NBS():" and words[2] == "WAVE[2Boct](start):"):
        tmp_time_sec[enum_nbs_oo] = conv_hmc_to_sec(words[9].strip())
    if (words[0] == "NBS_2Boct::run_NBS():" and words[2] == "WAVE[2Boct](end):"):
        num_count[enum_nbs_oo] += 1
        add_count(enum_nbs_oo, conv_hmc_to_sec(words[9].strip()))
    
    if (words[0] == "NBS_2Bdec::run_NBS():" and words[2] == "WAVE[2Bdec](start):"):
        tmp_time_sec[enum_nbs_dd] = conv_hmc_to_sec(words[9].strip())
    if (words[0] == "NBS_2Bdec::run_NBS():" and words[2] == "WAVE[2Bdec](end):"):
        num_count[enum_nbs_dd] += 1
        add_count(enum_nbs_dd, conv_hmc_to_sec(words[9].strip()))
    
    if (words[0] == "NBS_2Bdow::run_NBS():" and words[2] == "WAVE[2Bdow](start):"):
        tmp_time_sec[enum_nbs_od] = conv_hmc_to_sec(words[9].strip())
    if (words[0] == "NBS_2Bdow::run_NBS():" and words[2] == "WAVE[2Bdow](end):"):
        num_count[enum_nbs_od] += 1
        add_count(enum_nbs_od, conv_hmc_to_sec(words[9].strip()))

def conv_hmc_to_sec(hms):
    hms_sp = hms.split(":")
    return int(hms_sp[0])*60*60 + int(hms_sp[1])*60 + int(hms_sp[2])

def add_count(ienum, second):
    if (tmp_time_sec[ienum] > second):
        itmp = (24*60*60 - tmp_time_sec[ienum] + second)
    else:
        itmp = (second - tmp_time_sec[ienum])
    
    tot_time_sec [ienum] += itmp
    tot_time_sec2[ienum] += itmp**2

def make_mean_err_hms(ienum):
    sec_mean  = tot_time_sec [ienum] / float(num_count[ienum])
    sec_mean2 = tot_time_sec2[ienum] / float(num_count[ienum])
    if (num_count[ienum] != 1):
        sec_err = math.sqrt(abs((sec_mean2 - sec_mean**2)/float(num_count[ienum]-1)))
    else:
        sec_err = 0
    
    return (convert_sec_hms(int(sec_mean)), convert_sec_hms(int(sec_err)))

def convert_sec_hms(second):
    diff_sec  = (second % 3600) % 60
    diff_min  = ((second-diff_sec)/60) % 60
    diff_hour = (second-diff_min*60-diff_sec) / 3600
        
    return ("%02d:%02d:%02d" % (diff_hour, diff_min, diff_sec))

def print_results():
    print("Elapsed time:")
    for i in range(Ncount):
        if (num_count[i] != 0):
            tmp_res = make_mean_err_hms(i)
            print("%s = %s +/- %s (x%d)" % (enum_name[i], tmp_res[0], tmp_res[1], num_count[i]))

def check_results():
    tmp_sec = 0
    for i in range(1, Ncount):
        tmp_sec += tot_time_sec[i]
    
    print("Check:")
    print("SUM TIME = %s" % convert_sec_hms(tmp_sec))
    print("JOB TIME = %s" % convert_sec_hms(tot_time_sec[enum_totjob]))


if (__name__ == "__main__"):
    import sys
    ARGV = sys.argv; ARGC = len(ARGV)
    
    if (ARGC == 1):
        print("usage: %s [ifile(s)]" % ARGV[0]); quit()
    
    print("#.ifile = %d\n" % (ARGC-1))
    main(ARGC, ARGV)
    print_results()
    print
    check_results()
