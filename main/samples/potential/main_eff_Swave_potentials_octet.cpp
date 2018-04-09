#include <ComplexField_Sub.h>
#include <Potential.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int    Tsize = 64;
   int    Lsize = 32;
   double lattice_spacing = 0.0907;
   
   int bin_size = 57;
   
   int min_time_pot = 10;
   int max_time_pot = 12;
   
   bool input_from_bin_flg = true;
   bool output_corr_flg    = false;
   bool output_pot_flg     = true;
   
   string quark = "c";
   
   string ibase      = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/results_"+quark;
   string obase      = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/analysis/"+quark;
   string fconf_list = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/conf.lst";
   
   string ibase_corr = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/results.split/"+quark+"/results.split.ave";
   
   string odir_bin   = "bin";
   string odir_corr  = "corr";
   string odir_pot   = "pot1";
   
   string path_wave(string dir_base, string dir_channel, string conf, int time, string shift, string footer) {
      char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "%+04d", time);
      return (dir_base +"/"+ dir_channel +"/"+ conf +"/NBSwave."+ tmp_c +"+"+ shift +"."+ conf +"."+ footer);
   }
   string path_corr(string dir_base, string conf, string hadron, string shift) {
      return (dir_base +"/correlator.PS.dir/"+ conf + "/" + hadron +"_correlator.+"+ shift +"."+ conf);
   }
   
   string ichdir_str = "NBS_2Boct.Prot__Lamb__Prot__Lamb_.dir";
   
   string ihad_str[2] = {"proton_CG05_CG05", "Lambda_CG05_CG05"};
   
   string ch_str = "L"+quark+"N12__L"+quark+"N12_";
   
   string src_shift_str  = "Ave.000.000.000";
   string NBS_footer_str = "oo_CG05.oo_CG05";
   
   double mass_nucl  = 0.726840;
   double mass_Lam   = 0.754507;
   double mass_Lamc  = 1.233977;
   
   double tmp_mas_1, tmp_mas_2;
   
   anaHAL::NameList conf;
   STATISTICS<ComplexField_T  > corr[2];
   STATISTICS<ComplexField_XYZ> wave_1S0[3]; // 3 = t-1,t,t+1
   STATISTICS<ComplexField_XYZ> wave_3S1[3]; // 3 = t-1,t,t+1
   STATISTICS<ComplexField_XYZ> pot_eff_1S0, pot_eff_3S1;
   
   bool set_args(int, char**);
   
   void  input_corrdata_from_results();
   void  input_corrdata_from_bin    ();
   void output_corrdata_to_bin      ();
   void output_corrdata_to_text     ();
   
   void  input_NBSdata_from_results (const int);
   void  input_NBSdata_from_bin     (const int);
   void output_NBSdata_to_bin       (const int);
   
   void potential_calculation       (const int);
   void output_XYZdata_to_text      (const int);
   void output_XYZdata_to_bin       (const int);
   
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main(int argc, char **argv) {
   if(!set_args(argc, argv)) return -1;
   
   time_t stime, etime; time(&stime);
   
   if (input_from_bin_flg) input_corrdata_from_bin();
   else {
      input_corrdata_from_results();
      output_corrdata_to_bin();
   }
   
   if (output_corr_flg) output_corrdata_to_text();
   
   for (int itime=min_time_pot; itime<=max_time_pot; itime++) {
      printf("@@@ Potential calculation at t=%d START\n", itime);
      
      if (input_from_bin_flg) input_NBSdata_from_bin(itime);
      else {
         input_NBSdata_from_results(itime);
         output_NBSdata_to_bin(itime);
      }
      
      potential_calculation(itime);
      
      if (output_pot_flg) {
         output_XYZdata_to_text(itime);
         output_XYZdata_to_bin (itime);
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args(int argc_in, char **argv_in) {
      if      (quark == "c") {
         tmp_mas_1 = mass_nucl;
         tmp_mas_2 = mass_Lamc;
      }
      else if (quark == "s") {
         tmp_mas_1 = mass_nucl;
         tmp_mas_2 = mass_Lam ;
      }
      else {
         printf("ERROR quark, %s", quark.c_str());
         return false;
      }
      
      if      (argc_in == 1) return true;
      else if (argc_in != 7) {
         printf("usage: ./%s [ichdir_str] [ihad_str[0]] [ihad_str[1]] "
                "[mas_1] [mas_2] [ch_str(for output file name)]\n", argv_in[0]);
         return false;
      }
      else {
         ichdir_str  =      argv_in[1] ;
         ihad_str[0] =      argv_in[2] ;
         ihad_str[1] =      argv_in[3] ;
         tmp_mas_1   = atof(argv_in[4]);
         tmp_mas_2   = atof(argv_in[5]);
         ch_str      =      argv_in[6] ;
         return true;
      }
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_results() {
      conf.set(fconf_list);
      for (int ihad=0; ihad<2; ihad++) corr[ihad].mem_alloc(conf.Nlist());
      
      printf("Reading Correlators...\n");
      for (int ihad=0; ihad<2; ihad++) for (int i=0; i<conf.Nlist(); i++) {
         corr[ihad](i).mem_alloc(Tsize);
         corr[ihad](i).input_data_corr(path_corr(ibase_corr, conf(i), ihad_str[ihad], src_shift_str), true);
      }
      printf("Reading Correlators... end\n");
      
      printf("Make Jack-knife samples... \n");
      for (int ihad=0; ihad<2; ihad++) corr[ihad].make_JK_sample(bin_size);
      printf("Make Jack-knife samples... end\n");
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_bin() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      printf("Reading Correlators...\n");
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      printf("Reading Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_corrdata_to_bin() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      printf("Writing Correlators...\n");
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].output_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      printf("Writing Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_corrdata_to_text() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      printf("Writing Correlators...\n");
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].output_data_err(obase+"/"+odir_corr+"/correlator."+ihad_str[ihad]+
                                    ".bin_size"+Bsize_str+".gnu", 0, true);
      printf("Writing Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void input_NBSdata_from_results(const int itime) {
      conf.set(fconf_list);
      for (int dt=-1; dt<=+1; dt++) {
         wave_1S0[dt+1].mem_alloc(conf.Nlist());
         wave_3S1[dt+1].mem_alloc(conf.Nlist());
      }
      ComplexField_AXYZB tmp_wave(2*2, Lsize, 2*2);
      printf("Reading NBSwave...\n");
      for (int i=0; i<conf.Nlist(); i++) {
         for (int dt=-1; dt<=+1; dt++) {
            tmp_wave.input_data_bin(path_wave(ibase, ichdir_str, conf(i),
                                              itime+dt, src_shift_str, NBS_footer_str));
            wave_1S0[dt+1](i)  = tmp_wave.spin_proj(HH_OctOct, 0, 0, HH_OctOct, 0, 0);
            wave_3S1[dt+1](i)  = tmp_wave.spin_proj(HH_OctOct, 1,+1, HH_OctOct, 1,+1);
            wave_3S1[dt+1](i) += tmp_wave.spin_proj(HH_OctOct, 1, 0, HH_OctOct, 1, 0);
            wave_3S1[dt+1](i) += tmp_wave.spin_proj(HH_OctOct, 1,-1, HH_OctOct, 1,-1);
            wave_3S1[dt+1](i) /= 3.0;
         }
         printf("Reading NBSwave... end: conf=%d\n", i);
      }
      printf("Make Jack-knife samples... \n");
      for (int dt=-1; dt<=+1; dt++) {
         wave_1S0[dt+1].make_JK_sample(bin_size);
         wave_3S1[dt+1].make_JK_sample(bin_size);
      }
      printf("Make Jack-knife samples... end\n");
   }
   ////// ------------------------------------------------------ //////
   void input_NBSdata_from_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Reading NBSwave...\n");
      for (int dt=-1; dt<=+1; dt++) {
         char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
         string time_dt_str(time_dt_c);
         wave_1S0[dt+1].input_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_str+".spin_0_.t"+
                                       time_dt_str+".bin_size"+Bsize_str);
         wave_3S1[dt+1].input_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_str+".spin_1_.t"+
                                       time_dt_str+".bin_size"+Bsize_str);
      }
      printf("Reading NBSwave... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_NBSdata_to_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Writing NBSwave...\n");
      for (int dt=-1; dt<=+1; dt++) {
         char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
         string time_dt_str(time_dt_c);
         wave_1S0[dt+1].output_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_str+".spin_0_.t"+
                                        time_dt_str+".bin_size"+Bsize_str);
         wave_3S1[dt+1].output_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_str+".spin_1_.t"+
                                        time_dt_str+".bin_size"+Bsize_str);
      }
      printf("Writing NBSwave... end\n");
   }
   ////// ------------------------------------------------------ //////
   void potential_calculation(const int itime) {
      pot_eff_1S0.mem_alloc(wave_1S0[0].Ndata());
      pot_eff_3S1.mem_alloc(wave_3S1[0].Ndata());
      
      ComplexField_XYZ tmp_Rcor[3];
      printf("Calculating potentials... \n");
      for (int i=0; i<wave_1S0[0].Ndata(); i++) {
         for (int dt=-1; dt<=+1; dt++) tmp_Rcor[dt+1] = (wave_1S0[dt+1](i).rot_proj(ROT_REP_A1) /
                                                         (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
         
         pot_eff_1S0(i) = Potential::get_potential_T2(tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                                      tmp_mas_1, tmp_mas_2);
         
         for (int dt=-1; dt<=+1; dt++) tmp_Rcor[dt+1] = (wave_3S1[dt+1](i).rot_proj(ROT_REP_A1) /
                                                         (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
         
         pot_eff_3S1(i) = Potential::get_potential_T2(tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                                      tmp_mas_1, tmp_mas_2);
         
         printf("Calculating potentials... end: conf=%d\n", i);
      }
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_text(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      printf("Output potentials... \n");
      pot_eff_1S0.output_data_err(obase+"/"+odir_pot+"/Pot."+ch_str+".Veff1S0.t"+time_str+
                                  ".bin_size"+Bsize_str+".gnu", lattice_spacing, true);
      pot_eff_3S1.output_data_err(obase+"/"+odir_pot+"/Pot."+ch_str+".Veff3S1.t"+time_str+
                                  ".bin_size"+Bsize_str+".gnu", lattice_spacing, true);
      printf("Output potentials... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      printf("Output potentials... \n");
      pot_eff_1S0.output_data_bin_reduce(obase+"/"+odir_pot+"/Pot."+ch_str+".Veff1S0.t"+time_str+
                                         ".bin_size"+Bsize_str+".bin", lattice_spacing, false);
      pot_eff_3S1.output_data_bin_reduce(obase+"/"+odir_pot+"/Pot."+ch_str+".Veff3S1.t"+time_str+
                                         ".bin_size"+Bsize_str+".bin", lattice_spacing, false);
      printf("Output potentials... end\n");
   }
}
