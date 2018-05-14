#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>

namespace  {
   int    Tsize = 32;
   int    Lsize = 16;
   double lattice_spacing = 0.1215;
   
   int bin_size = 35;
   
   string ibase_rl   = "/xc/home/takaya.miyamoto/work/su3_16x32/run_rela/results";
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
   string tmp_str    [4] = {"NN"   , "ND"   , "DN"   , "DD"   };
   string tmp_str_org[4] = {"NN_rl", "NN_nr", "DD_rl", "DD_nr"};
   
   string ch_decdec[2] = {"DelPP_DelM_", "DelP__DelZ_"};
   
   string sft  = "Ave.000.000.000";
   
   double mass_nuc = 1.242;
   double mass_del = 1.363;
   
   double ZP_nuc_rl = pow(0.316385, 2);
   double ZP_nuc_nr = pow(0.246395, 2);
   double ZP_del_rl = pow(0.388392, 2);
   double ZP_del_nr = pow(0.0111584, 2) * 1024; // Adjustment factor for NR-Delta <--> NR-OMEGA (Yamada-san's definition)
}

int main() {
   time_t stime, etime; time(&stime);
   
   //anaHAL::NameList conf(fconf_list);
   //printf("#.conf = %d\n", conf.Nlist());
   
   STATISTICS<ComplexField_T  > corr_NN_nr;//(conf.Nlist());
   STATISTICS<ComplexField_T  > corr_NN_rl;//(conf.Nlist());
   STATISTICS<ComplexField_T  > corr_DD_nr;//(conf.Nlist());
   STATISTICS<ComplexField_T  > corr_DD_rl;//(conf.Nlist());
   STATISTICS<ComplexField_XYZ> wave_NN_nr[3];
   STATISTICS<ComplexField_XYZ> wave_NN_rl[3];
   STATISTICS<ComplexField_XYZ> wave_DD_nr[3];
   STATISTICS<ComplexField_XYZ> wave_DD_rl[3];
   
   ComplexField_AXYZB wave_ooxdd(2*2,     Lsize, 2*2*3*3);
   ComplexField_AXYZB wave_ddxdd(2*2*3*3, Lsize, 2*2*3*3);
   ComplexField_XYZ   tmp_wave  [4];
   ComplexField_XYZ Y20      = sfunc::cfield_Ylm(2, 0, Lsize);
   ComplexField_XYZ Y20_2_A1 = Y20*Y20; Y20_2_A1.rot_proj(ROT_REP_A1);
   
   char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "%d", bin_size); string Bsize_str(tmp_c);
   /*
    printf("Reading Correlators...\n");
    for (int i=0; i<corr_NN_nr.Ndata(); i++) {
    corr_NN_nr(i).mem_alloc(Tsize); corr_NN_rl(i).mem_alloc(Tsize);
    corr_DD_nr(i).mem_alloc(Tsize); corr_DD_rl(i).mem_alloc(Tsize);
    corr_NN_nr(i).input_data_corr(path_corr(ibase_nr, conf(i), "proton_NR_NR"    , sft), true);
    corr_NN_rl(i).input_data_corr(path_corr(ibase_nr, conf(i), "proton_CG05_CG05", sft), true);
    corr_DD_nr(i).input_data_corr(path_corr(ibase_nr, conf(i), "OMEGA"           , sft), true);
    corr_DD_rl(i).input_data_corr(path_corr(ibase_nr, conf(i), "Delta"           , sft), true);
    }
    printf("Reading Correlators... end\n");
    
    printf("Make Jack-knife samples... \n");
    corr_NN_rl.make_JK_sample(bin_size);
    corr_NN_nr.make_JK_sample(bin_size);
    corr_DD_rl.make_JK_sample(bin_size);
    corr_DD_nr.make_JK_sample(bin_size);
    corr_DD_nr *= 1024; // Adjustment factor for NR-Delta <--> NR-OMEGA (Yamada-san's definition)
    printf("Make Jack-knife samples... end\n");
    
    corr_NN_rl.output_data_bin(obase+"/bin/correlator"+"."+tmp_str_org[0]+".bin_size"+Bsize_str);
    corr_NN_nr.output_data_bin(obase+"/bin/correlator"+"."+tmp_str_org[1]+".bin_size"+Bsize_str);
    corr_DD_rl.output_data_bin(obase+"/bin/correlator"+"."+tmp_str_org[2]+".bin_size"+Bsize_str);
    corr_DD_nr.output_data_bin(obase+"/bin/correlator"+"."+tmp_str_org[3]+".bin_size"+Bsize_str);
    */
   corr_NN_rl.input_data_bin(ibase+"/bin/correlator"+"."+tmp_str_org[0]+".bin_size"+Bsize_str);
   corr_NN_nr.input_data_bin(ibase+"/bin/correlator"+"."+tmp_str_org[1]+".bin_size"+Bsize_str);
   corr_DD_rl.input_data_bin(ibase+"/bin/correlator"+"."+tmp_str_org[2]+".bin_size"+Bsize_str);
   corr_DD_nr.input_data_bin(ibase+"/bin/correlator"+"."+tmp_str_org[3]+".bin_size"+Bsize_str);
   
   for (int itime=7; itime<11; itime++) {
      printf("@@@ NBSwave t=%d START\n", itime);
      char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "t%d", itime); string time_str(tmp_c);
      /*
       for (int i=0; i<3; i++) {
       wave_NN_nr[i].mem_alloc(conf.Nlist());
       wave_NN_rl[i].mem_alloc(conf.Nlist());
       wave_DD_nr[i].mem_alloc(conf.Nlist());
       wave_DD_rl[i].mem_alloc(conf.Nlist());
       }
       
       printf("Reading NBSwave...\n");
       for (int i=0; i<wave_NN_rl[0].Ndata(); i++) {
       for (int dt=-1; dt<=+1; dt++) {
       
       for (int j=0; j<2; j++) for (int k=0; k<2; k++) {
       wave_ddxdd.input_data_bin(path_wave(ibase_rl, ("NBS_2Bdec."+ch_decdec[j]+"_"+ch_decdec[k]+".dir"),
       conf(i), itime+dt, sft, "dd_CG0K.dd_CG0K"));
       tmp_wave[k+2*j] = wave_ddxdd.spin_proj(HH_DecDec, 3, 0, HH_DecDec, 3, 0);
       }
       wave_DD_rl[dt+1](i) = (tmp_wave[0] - tmp_wave[1] - tmp_wave[2] + tmp_wave[3]) / 4.0;
       
       for (int j=0; j<2; j++) for (int k=0; k<2; k++)
       tmp_wave[k+2*j].input_data_bin(path_wave(ibase_nr,
       ("Proj.NBS_2Bdec."+ch_decdec[j]+"_"+ch_decdec[k]+".dir/spin3z+0.3z+0"),
       conf(i), itime+dt, sft, "dd_NR.dd_NR"));
       wave_DD_nr[dt+1](i) = (tmp_wave[0] - tmp_wave[1] - tmp_wave[2] + tmp_wave[3]) / 4.0;
       
       for (int k=0; k<2; k++) {
       wave_ooxdd.input_data_bin(path_wave(ibase_rl, ("NBS_ooxdd.Prot__Neut__"+ch_decdec[k]+".dir"),
       conf(i), itime+dt, sft, "oo_CG05.dd_CG0K"));
       tmp_wave[k] = wave_ooxdd.spin_proj(HH_OctOct, 1, 0, HH_DecDec, 3, 0);
       }
       wave_NN_rl[dt+1](i) = (tmp_wave[0] - tmp_wave[1]) / 2.0;
       
       for (int k=0; k<2; k++)
       tmp_wave[k].input_data_bin(path_wave(ibase_nr,
       ("Proj.NBS_ooxdd.Prot__Neut__"+ch_decdec[k]+".dir/spin1z+0.3z+0"),
       conf(i), itime+dt, sft, "oo_NR.dd_NR"));
       wave_NN_nr[dt+1](i) = (tmp_wave[0] - tmp_wave[1]) / 2.0;
       }
       printf("Reading NBSwave... end: conf=%d\n", i);
       
       } // iconf
       */
      printf("Make Jack-knife samples... \n");
      for (int dt=-1; dt<=+1; dt++) {
         char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "t%d", itime+dt); string time_str_dt(tmp_c);
         /*
          wave_NN_rl[dt+1].make_JK_sample(bin_size);
          wave_NN_nr[dt+1].make_JK_sample(bin_size);
          wave_DD_rl[dt+1].make_JK_sample(bin_size);
          wave_DD_nr[dt+1].make_JK_sample(bin_size);
          
          wave_NN_rl[dt+1].output_data_bin(obase+"/bin/NBSwave"+"."+tmp_str_org[0]+"."+time_str_dt+".bin_size"+Bsize_str);
          wave_NN_nr[dt+1].output_data_bin(obase+"/bin/NBSwave"+"."+tmp_str_org[1]+"."+time_str_dt+".bin_size"+Bsize_str);
          wave_DD_rl[dt+1].output_data_bin(obase+"/bin/NBSwave"+"."+tmp_str_org[2]+"."+time_str_dt+".bin_size"+Bsize_str);
          wave_DD_nr[dt+1].output_data_bin(obase+"/bin/NBSwave"+"."+tmp_str_org[3]+"."+time_str_dt+".bin_size"+Bsize_str);
          */
         wave_NN_rl[dt+1].input_data_bin(ibase+"/bin/NBSwave"+"."+tmp_str_org[0]+"."+time_str_dt+".bin_size"+Bsize_str);
         wave_NN_nr[dt+1].input_data_bin(ibase+"/bin/NBSwave"+"."+tmp_str_org[1]+"."+time_str_dt+".bin_size"+Bsize_str);
         wave_DD_rl[dt+1].input_data_bin(ibase+"/bin/NBSwave"+"."+tmp_str_org[2]+"."+time_str_dt+".bin_size"+Bsize_str);
         wave_DD_nr[dt+1].input_data_bin(ibase+"/bin/NBSwave"+"."+tmp_str_org[3]+"."+time_str_dt+".bin_size"+Bsize_str);
      }
      printf("Make Jack-knife samples... end\n");
      
      STATISTICS<ComplexField_XYZ> detR(wave_NN_rl[0].Ndata());
      
      STATISTICS<ComplexField_XYZ> SKR[4];
      STATISTICS<ComplexField_XYZ> pot[4];
      for (int k=0; k<4; k++) {
         SKR[k].mem_alloc(wave_NN_rl[0].Ndata());
         pot[k].mem_alloc(wave_NN_rl[0].Ndata());
      }
      
      ComplexFieldMatrix<ComplexField_XYZ>  R(2);
      ComplexFieldMatrix<ComplexField_XYZ> KR(2);
      ComplexFieldMatrix<ComplexField_XYZ>  P(2);
      
      printf("Calculating potentials... \n");
      for (int i=0; i<wave_NN_rl[0].Ndata(); i++) {
         for (int dt=-1; dt<=+1; dt++) {
            wave_NN_nr[dt+1](i) /= (corr_NN_nr(i)(itime+dt) * corr_NN_nr(i)(itime+dt));
            wave_NN_rl[dt+1](i) /= (corr_NN_rl(i)(itime+dt) * corr_NN_rl(i)(itime+dt));
            wave_DD_nr[dt+1](i) /= (corr_DD_nr(i)(itime+dt) * corr_DD_nr(i)(itime+dt));
            wave_DD_rl[dt+1](i) /= (corr_DD_rl(i)(itime+dt) * corr_DD_rl(i)(itime+dt));
            
            wave_DD_nr[dt+1](i)  = wave_DD_nr[dt+1](i).rot_proj(ROT_REP_A1);
            wave_DD_rl[dt+1](i)  = wave_DD_rl[dt+1](i).rot_proj(ROT_REP_A1);
            
            //wave_NN_nr[dt+1](i)  = wave_NN_nr[dt+1](i) - wave_NN_nr[dt+1](i).rot_proj(ROT_REP_A1);
            //wave_NN_rl[dt+1](i)  = wave_NN_nr[dt+1](i) - wave_NN_rl[dt+1](i).rot_proj(ROT_REP_A1);
            
            wave_NN_nr[dt+1](i)  = wave_NN_nr[dt+1](i).rot_proj(ROT_REP_E);
            wave_NN_rl[dt+1](i)  = wave_NN_rl[dt+1](i).rot_proj(ROT_REP_E);
            wave_NN_nr[dt+1](i) *= Y20; wave_NN_nr[dt+1](i).rot_proj(ROT_REP_A1); wave_NN_nr[dt+1](i) /= Y20_2_A1;
            wave_NN_rl[dt+1](i) *= Y20; wave_NN_rl[dt+1](i).rot_proj(ROT_REP_A1); wave_NN_rl[dt+1](i) /= Y20_2_A1;
         }
         cdouble del_nr = ((ZP_nuc_nr * corr_DD_nr(i)(itime) * corr_DD_nr(i)(itime)) /
                           (ZP_del_nr * corr_NN_nr(i)(itime) * corr_NN_nr(i)(itime)));
         cdouble del_rl = ((ZP_nuc_rl * corr_DD_rl(i)(itime) * corr_DD_rl(i)(itime)) /
                           (ZP_del_rl * corr_NN_rl(i)(itime) * corr_NN_rl(i)(itime)));
         
         KR(0,0) = (wave_NN_rl[1](i).lap() / mass_nuc + (wave_NN_rl[0](i) - wave_NN_rl[2](i)) / 2.0 +
                    (wave_NN_rl[0](i) + wave_NN_rl[2](i) - 2.0 * wave_NN_rl[1](i)) / (4.0 * mass_nuc));
         
         KR(0,1) = (wave_NN_nr[1](i).lap() / mass_nuc + (wave_NN_nr[0](i) - wave_NN_nr[2](i)) / 2.0 +
                    (wave_NN_nr[0](i) + wave_NN_nr[2](i) - 2.0 * wave_NN_nr[1](i)) / (4.0 * mass_nuc));
         
         KR(1,0) = (wave_DD_rl[1](i).lap() / mass_del + (wave_DD_rl[0](i) - wave_DD_rl[2](i)) / 2.0 +
                    (wave_DD_rl[0](i) + wave_DD_rl[2](i) - 2.0 * wave_DD_rl[1](i)) / (4.0 * mass_del)) * del_rl;
         
         KR(1,1) = (wave_DD_nr[1](i).lap() / mass_del + (wave_DD_nr[0](i) - wave_DD_nr[2](i)) / 2.0 +
                    (wave_DD_nr[0](i) + wave_DD_nr[2](i) - 2.0 * wave_DD_nr[1](i)) / (4.0 * mass_del)) * del_nr;
         
         wave_DD_rl[1](i) *= del_rl;
         wave_DD_nr[1](i) *= del_nr;
         
         R (0,0) = wave_NN_rl[1](i);
         R (0,1) = wave_NN_nr[1](i);
         R (1,0) = wave_DD_rl[1](i);
         R (1,1) = wave_DD_nr[1](i);
         
         P       = KR * R.inverce();
         detR(i) = R.det();
         
         for (int k=0; k<4; k++) {
            SKR[k](i) = KR(k);
            pot[k](i) =  P(k);
         }
         printf("Calculating potentials... end: conf=%d\n", i);
         
      } // iconf
      
      string ofname_base;
      printf("Output potentials... \n");
      
      ofname_base = "/wave/Rcorr";
      
      wave_NN_rl[1].output_data_err(obase+ofname_base+"."+tmp_str[0]+"."+time_str+".gnu", 0, true);
      wave_NN_nr[1].output_data_err(obase+ofname_base+"."+tmp_str[1]+"."+time_str+".gnu", 0, true);
      wave_DD_rl[1].output_data_err(obase+ofname_base+"."+tmp_str[2]+"."+time_str+".gnu", 0, true);
      wave_DD_nr[1].output_data_err(obase+ofname_base+"."+tmp_str[3]+"."+time_str+".gnu", 0, true);
      
      ofname_base = "/wave/Rcorr";
      detR.output_data_err(obase+ofname_base+".det"+"."+time_str+".gnu", 0, true);
      
      ofname_base = "/pot/Pot";
      for (int i=0; i<4; i++)
         pot[i].output_data_err(obase+ofname_base+"."+tmp_str[i]+"."+time_str+".gnu", lattice_spacing, true);
      
      ofname_base = "/pot/KerR";
      for (int i=0; i<4; i++)
         SKR[i].output_data_err(obase+ofname_base+"."+tmp_str[i]+"."+time_str+".gnu", 0, true);
      
      printf("Output potentials... end\n");
      printf("@@@ NBSwave t=%d  END \n", itime);
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
