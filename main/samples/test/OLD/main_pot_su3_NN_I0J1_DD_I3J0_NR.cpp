#include <ComplexField_Sub.h>

namespace  {
   int    Tsize = 32;
   int    Lsize = 16;
   double lattice_spacing = 0.1215;
   
   int bin_size = 35;
   
   string ibase_nr   = "/xc/home/takaya.miyamoto/data/su3/16x32/results_NR";
   string ibase      = "/home/soryushi/takaya.miyamoto/xchome/work/su3_16x32/analysis";
   string obase      = "/home/soryushi/takaya.miyamoto/xchome/work/su3_16x32/analysis_new";
   string fconf_list = obase + "/conf.lst";
   
   string path_wave(string dir_base, string dir_channel, string conf, int time, string shift, string footer) {
      char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "%+04d", time);
      return (dir_base +"/"+ dir_channel +"/"+ conf +"/NBSwave."+ tmp_c +"+"+ shift +"."+ conf +"."+ footer);
   }
   string path_corr(string dir_base, string conf, string hadron, string shift) {
      return (dir_base +"/correlator.PS.dir/"+ conf + "/" + hadron +"_correlator.+"+ shift +"."+ conf);
   }
   
   string ch_decdec[2] = {"DelPP_DelM_", "DelP__DelZ_"};
   
   string tmp_str_org[4] = {"NN_rl", "NN_nr", "DD_rl", "DD_nr"};
   
   string sft = "Ave.000.000.000";
   
   double mass_nuc = 1.242;
   double mass_del = 1.363;
   
   string ofname_base;
}

int main() {
   time_t stime, etime; time(&stime);
   
   //anaHAL::NameList conf(fconf_list);
   //printf("#.conf = %d\n", conf.Nlist());
   
   STATISTICS<ComplexField_T  > corr_NN_nr;//(conf.Nlist());
   STATISTICS<ComplexField_T  > corr_DD_nr;//(conf.Nlist());
   
   STATISTICS<ComplexField_XYZ> wave_NN_nr[3];
   STATISTICS<ComplexField_XYZ> wave_DD_nr[3];
   
   ComplexField_XYZ Y20      = sfunc::cfield_Ylm(2, 0, Lsize);
   ComplexField_XYZ Y20_2_A1 = Y20*Y20; Y20_2_A1.rot_proj(ROT_REP_A1);
   
   ComplexField_XYZ tmp_wave[4];
   
   for (int i=0; i<4; i++) tmp_wave[i].mem_alloc(Lsize);
   
   char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "%d", bin_size); string Bsize_str(tmp_c);
   /*
    printf("Reading Correlators...\n");
    for (int i=0; i<corr_NN_nr.Ndata(); i++) {
    corr_NN_nr(i).mem_alloc(Tsize);
    corr_DD_nr(i).mem_alloc(Tsize);
    corr_NN_nr(i).input_data_corr(path_corr(ibase_nr, conf(i), "proton_NR_NR", sft), true);
    corr_DD_nr(i).input_data_corr(path_corr(ibase_nr, conf(i), "OMEGA"       , sft), true);
    }
    printf("Reading Correlators... end\n");
    
    printf("Make Jack-knife samples... \n");
    corr_NN_nr.make_JK_sample(bin_size);
    corr_DD_nr.make_JK_sample(bin_size);
    corr_DD_nr *= 1024; // Adjustment factor for NR-Delta <--> NR-OMEGA (Yamada-san's definition)
    printf("Make Jack-knife samples... end\n");
    */
   //corr_NN_nr.output_data_bin(obase+"/bin/correlator"+"."+tmp_str_org[1]+".bin_size"+Bsize_str);
   //corr_DD_nr.output_data_bin(obase+"/bin/correlator"+"."+tmp_str_org[3]+".bin_size"+Bsize_str);
   corr_NN_nr.input_data_bin(ibase+"/bin/correlator"+"."+tmp_str_org[1]+".bin_size"+Bsize_str);
   corr_DD_nr.input_data_bin(ibase+"/bin/correlator"+"."+tmp_str_org[3]+".bin_size"+Bsize_str);
   
   ofname_base = "/corr/corr";
   corr_NN_nr.output_data_err(obase+ofname_base+".NN_nr.gnu", 0, true);
   corr_DD_nr.output_data_err(obase+ofname_base+".DD_nr.gnu", 0, true);
   
   for (int itime=7; itime<11; itime++) {
      printf("@@@ NBSwave t=%d START\n", itime);
      char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "t%d", itime); string time_str(tmp_c);
      /*
       for (int i=0; i<3; i++) {
       wave_NN_nr[i].mem_alloc(conf.Nlist());
       wave_DD_nr[i].mem_alloc(conf.Nlist());
       }
       
       printf("Reading NBSwave...\n");
       for (int i=0; i<wave_NN_nr[0].Ndata(); i++) {
       for (int dt=-1; dt<=+1; dt++) {
       
       for (int j=0; j<2; j++) for (int k=0; k<2; k++)
       tmp_wave[k+2*j].input_data_bin(path_wave(ibase_nr,
       ("Proj.NBS_2Bdec."+ch_decdec[j]+"_"+ch_decdec[k]+".dir/spin3z+0.3z+0"),
       conf(i), itime+dt, sft, "dd_NR.dd_NR"));
       wave_DD_nr[dt+1](i) = (tmp_wave[0] - tmp_wave[1] - tmp_wave[2] + tmp_wave[3]) / 4.0;
       
       tmp_wave[0].input_data_bin(path_wave(ibase_nr, "Proj.NBS_2Boct.Prot__Neut__Prot__Neut_.dir/spin1z+0.1z+0",
       conf(i), itime+dt, sft, "oo_NR.oo_NR"));
       wave_NN_nr[dt+1](i) = tmp_wave[0];
       }
       printf("Reading NBSwave... end: conf=%d\n", i);
       
       } // iconf
       */
      printf("Make Jack-knife samples... \n");
      for (int dt=-1; dt<=+1; dt++) {
         char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "t%d", itime+dt); string time_str_dt(tmp_c);
         /*
          wave_NN_nr[dt+1].make_JK_sample(bin_size);
          wave_DD_nr[dt+1].make_JK_sample(bin_size);
          wave_NN_nr[dt+1].output_data_bin(obase+"/bin/NBSwave.NN_nr_single."+time_str_dt+".bin_size"+Bsize_str);
          wave_DD_nr[dt+1].output_data_bin(obase+"/bin/NBSwave.DD_nr_single."+time_str_dt+".bin_size"+Bsize_str);
          */
         wave_NN_nr[dt+1].input_data_bin(ibase+"/bin/NBSwave.NN_nr_single."+time_str_dt+".bin_size"+Bsize_str);
         wave_DD_nr[dt+1].input_data_bin(ibase+"/bin/NBSwave.DD_nr_single."+time_str_dt+".bin_size"+Bsize_str);
      }
      printf("Make Jack-knife samples... end\n");
      
      STATISTICS<ComplexField_XYZ> pot_NN(wave_NN_nr[0].Ndata());
      STATISTICS<ComplexField_XYZ> pot_DD(wave_NN_nr[0].Ndata());
      
      printf("Calculating potentials... \n");
      for (int i=0; i<wave_NN_nr[0].Ndata(); i++) {
         for (int dt=-1; dt<=+1; dt++) {
            wave_NN_nr[dt+1](i) /= (corr_NN_nr(i)(itime+dt) * corr_NN_nr(i)(itime+dt));
            wave_DD_nr[dt+1](i) /= (corr_DD_nr(i)(itime+dt) * corr_DD_nr(i)(itime+dt));
            
            wave_DD_nr[dt+1](i)  = wave_DD_nr[dt+1](i).rot_proj(ROT_REP_A1);
            
            //wave_NN_nr[dt+1](i)  = wave_NN_nr[dt+1](i) - wave_NN_nr[dt+1](i).rot_proj(ROT_REP_A1);
            
            wave_NN_nr[dt+1](i)  = wave_NN_nr[dt+1](i).rot_proj(ROT_REP_E);
            wave_NN_nr[dt+1](i) *= Y20; wave_NN_nr[dt+1](i).rot_proj(ROT_REP_A1); wave_NN_nr[dt+1](i) /= Y20_2_A1;
         }
         pot_NN(i) = (wave_NN_nr[1](i).lap() / mass_nuc + (wave_NN_nr[0](i) - wave_NN_nr[2](i)) / 2.0 +
                      (wave_NN_nr[0](i) + wave_NN_nr[2](i) - 2.0 * wave_NN_nr[1](i)) / (4.0 * mass_nuc)) / wave_NN_nr[1](i);
         pot_DD(i) = (wave_DD_nr[1](i).lap() / mass_del + (wave_DD_nr[0](i) - wave_DD_nr[2](i)) / 2.0 +
                      (wave_DD_nr[0](i) + wave_DD_nr[2](i) - 2.0 * wave_DD_nr[1](i)) / (4.0 * mass_del)) / wave_DD_nr[1](i);
         
         printf("Calculating potentials... end: conf=%d\n", i);
         
      } // iconf
      
      printf("Output potentials... \n");
      
      ofname_base = "/wave/Rcorr";
      wave_NN_nr[1].output_data_err(obase+ofname_base+".NN_nr_single."+time_str+".gnu", 0, true);
      wave_DD_nr[1].output_data_err(obase+ofname_base+".DD_nr_single."+time_str+".gnu", 0, true);
      
      ofname_base = "/pot/Pot";
      pot_NN.output_data_err(obase+ofname_base+".NN_nr_single."+time_str+".gnu", lattice_spacing, true);
      pot_DD.output_data_err(obase+ofname_base+".DD_nr_single."+time_str+".gnu", lattice_spacing, true);
      
      printf("Output potentials... end\n");
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
