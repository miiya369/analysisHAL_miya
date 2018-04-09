#!/usr/bin/python

Tsize = 32
Ncomp = 11

ifbase = "/home/soryushi/takaya.miyamoto/xchome/work/cp-pacs.faisal/work.run/results.split"
ofbase = "/home/soryushi/takaya.miyamoto/data/cp-pacs.faisal"
ifclst = "/home/soryushi/takaya.miyamoto/xchome/work/cp-pacs.faisal/work.run/results/conf_all.lst"
ofclst = "/home/soryushi/takaya.miyamoto/data/cp-pacs.faisal/conf_out.lst"
fpsss = ifbase + "/correlator.ps_ss.lst"
fmult = ifbase + "/correlator.multi.lst"

if (__name__ == "__main__"):
    import os, sys
    import numpy as np
    
    if (len(sys.argv) != 2):
        print("usage: %s [time%%4 (=0,1,2,3)]" % sys.argv[0]); quit()
    t_mod_04 = int(sys.argv[1].strip()) #; print t_mod_04; quit()
    
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
    
    for iconf_out in range(Nconf_out):
        for corr_dir in ("correlator.PS.dir", "correlator.SS.dir", "correlator.multi.dir"):
            tmpTsize  = Tsize / 2 if (corr_dir == "correlator.multi.dir") else Tsize
            fhadname  = fmult     if (corr_dir == "correlator.multi.dir") else fpsss
            odir_name = ("%s/t_mod04_%03d/%s/%s" % 
                         (ofbase, t_mod_04, corr_dir, conf_lst_out[iconf_out].strip())); 
            os.makedirs(odir_name)
            
            for had_name in open(fhadname, 'r'):
                corr = np.zeros((tmpTsize, Ncomp))
                itmp = 0
                for iconf_in in range(iconf_out * Nconf_ave, (iconf_out + 1) * Nconf_ave):
                    for tsrc in range(t_mod_04, Tsize, 4):
                        ifile_name = ("%s/%s/%s/%s.+%03d.000.000.000.%s" % 
                                      (ifbase, corr_dir, conf_lst_in[iconf_in].strip(), 
                                       had_name.strip(), tsrc, conf_lst_in[iconf_in].strip()))
                        tmp_corr = np.loadtxt(ifile_name)
                        corr    += tmp_corr
                        itmp    += 1; del tmp_corr
                if (itmp != Naverage):
                    print("\nERROR: Unexpected Naverage (Total ave. != %d), exit\n" % itmp); quit()
                corr /= float(Naverage)
                
                ofile_name = ("%s/%s.+A%02d.000.000.000.%s" %
                              (odir_name, had_name.strip(), t_mod_04, conf_lst_out[iconf_out].strip()))
                np.savetxt(ofile_name, corr); del corr
                print(" - End: conf = %03d, dir = %s, had = %s" % (iconf_out, corr_dir, had_name.strip()))
