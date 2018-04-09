#!/usr/bin/python

Tsize = 32
Lsize = 16

t_min = 0
t_max = 15

ifbase = "/home/soryushi/takaya.miyamoto/xchome/work/cp-pacs.faisal/work.run/results"
ofbase = "/home/soryushi/takaya.miyamoto/data/cp-pacs.faisal"
ifclst = "/home/soryushi/takaya.miyamoto/xchome/work/cp-pacs.faisal/work.run/results/conf_all.lst"
ofclst = "/home/soryushi/takaya.miyamoto/data/cp-pacs.faisal/conf_out.lst"

nbs_footer = "op_CG05.op_CG05"

if (__name__ == "__main__"):
    import os, sys
    import numpy as np
    from struct import pack
    
    if (len(sys.argv) < 3):
        print("usage: %s [time%%4 (=0,1,2,3)] [NBS dir1] [NBS dir2] ..." % sys.argv[0]); quit()
    t_mod_04 = int(sys.argv[1].strip())
    Nnbs_dir = len(sys.argv) - 2; nbs_dir = [None] * Nnbs_dir
    for i in range(Nnbs_dir):
        nbs_dir[i] = sys.argv[i+2].strip()
    #print t_mod_04, nbs_dir; quit()
    
    with open(ifclst, 'r') as cf:
        conf_lst_in  = cf.readlines() #; print conf_lst_in [6]; quit()
        Nconf_in     = len(conf_lst_in)
    with open(ofclst, 'r') as cf:
        conf_lst_out = cf.readlines() #; print conf_lst_out[6]; quit()
        Nconf_out    = len(conf_lst_out)
    
    if (Nconf_in % Nconf_out != 0):
        print("\nERROR: Nconf_in = %d cannot be devided by Nconf_out = %d, exit.\n" % 
              (Nconf_in, Nconf_out)); quit()
    
    Nconf_ave = Nconf_in / Nconf_out
    Naverage  = Nconf_ave * 8
    print("\n# Nconf_in = %d, Nconf_out = %d, Nconf_ave = %d (Total ave. = %d): Average Start.\n" % 
          (Nconf_in, Nconf_out, Nconf_ave, Naverage)) #; quit()
    
    nbs_size = Lsize**3 * 2
    
    for iconf_out in range(Nconf_out):
        for inbs  in range(Nnbs_dir):
            odir_name = ("%s/t_mod04_%03d/%s/%s" % 
                         (ofbase, t_mod_04, nbs_dir[inbs], conf_lst_out[iconf_out].strip())); 
            os.makedirs(odir_name)
            
            for it in range(t_min, t_max+1):
                nbs  = np.zeros(nbs_size)
                itmp = 0
                for iconf_in in range(iconf_out * Nconf_ave, (iconf_out + 1) * Nconf_ave):
                    for tsrc in range(t_mod_04, Tsize, 4):
                        ifile_name = ("%s/%s/%s/NBSwave.+%03d+%03d.000.000.000.%s.%s" % 
                                      (ifbase, nbs_dir[inbs], conf_lst_in[iconf_in].strip(), it, tsrc, 
                                       conf_lst_in[iconf_in].strip(), nbs_footer))
                        tmp_nbs = np.fromfile(ifile_name, '>d') #; print tmp_nbs; quit()
                        nbs    += tmp_nbs
                        itmp   += 1; del tmp_nbs
                if (itmp != Naverage):
                    print("\nERROR: Unexpected Naverage (Total ave. != %d), exit\n" % itmp); quit()
                nbs /= float(Naverage)
                
                ofile_name = ("%s/NBSwave.+%03d+A%02d.000.000.000.%s.%s" %
                              (odir_name, it, t_mod_04, conf_lst_out[iconf_out].strip(), nbs_footer))
                with open(ofile_name, 'wb') as fnbs:
                    for n in range(nbs_size):
                        fnbs.write(pack('>d', nbs[n]))
                del nbs
                print(" - End: conf = %03d, dir = %s, time = %03d" % (iconf_out, nbs_dir[inbs], it))
