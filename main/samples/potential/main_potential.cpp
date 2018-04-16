#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>
#include <Potential.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int    Tsize = 64;
   double lattice_spacing = 0.0907;
   
   int bin_size = 57;
   
   int min_time_pot = 10;
   int max_time_pot = 12;
   
   bool pre_read_flg = true;
   
   string quark = "c";
   
   string ibase      = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/results_"+quark;
   string obase      = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/analysis/"+quark;
   string fconf_list = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/conf.lst";
   
   string ibase_corr = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela/results.split/"+quark+"/results.split.ave";
   
   string odir_bin   = "bin";
   string odir_corr  = "corr";
   string odir_pot   = "pot3";
   
   string path_corr(string dir_base, string conf, string hadron, string shift) {
      return (dir_base +"/correlator.PS.dir/"+ conf + "/" + hadron +"_correlator.+"+ shift +"."+ conf);
   }
   
   string ihad_3x3_str[3][2] = {
      {"proton_CG05_CG05", "Lambda_CG05_CG05"},
      {"proton_CG05_CG05", "Sigma_CG05_CG05"},
      {"proton_CG05_CG05", "Sigma32"}
   };
   string ihad_2x2_str[2][2] = {
      {"proton_CG05_CG05", "Sigma_CG05_CG05"},
      {"proton_CG05_CG05", "Sigma32"}
   };
   
   string ch_3x3_str[3] = {"L"+quark+"N12_", "S"+quark+"N12_", "S"+quark+"SN12"};
   string ch_2x2_str[2] = {"S"+quark+"N32_", "S"+quark+"SN32"};
  
   string src_shift_str  = "A64.000.000.000";
   
   double mass_nucl  = 0.726840;
   double mass_Lam   = 0.754507;
   double mass_Sig   = 0.761660;
   double mass_SigS  = 0.861298;
   double mass_Lamc  = 1.233977;
   double mass_Sigc  = 1.277982;
   double mass_SigcS = 1.317399;
   
   double sq_ZP_nucl  = 0.0638732520686627;
   double sq_ZP_Lam   = 0.0681264891598947;
   double sq_ZP_Sig   = 0.0672251624030046;
   double sq_ZP_SigS  = 0.0784868204517082;
   double sq_ZP_Lamc  = 0.1176390522262307;
   double sq_ZP_Sigc  = 0.1025945201743403;
   double sq_ZP_SigcS = 0.1402728573790395;
   
   double tmp_mas_3x3_1[3], tmp_mas_3x3_2[3];
   double tmp_mas_2x2_1[2], tmp_mas_2x2_2[2];
   
   cdouble ZP_3x3[3], ZP_2x2[2];
   
   anaHAL::NameList conf;
   STATISTICS<ComplexField_T  > corr_3x3[3][2];
   STATISTICS<ComplexField_T  > corr_2x2[2][2];
   
   STATISTICS<ComplexField_XYZ> wave_3S1_3x3[3][3*3];
   STATISTICS<ComplexField_XYZ> wave_3S1_2x2[3][2*2];
   
   STATISTICS<ComplexField_XYZ> pot_eff_3S1_3x3[3*3];
   STATISTICS<ComplexField_XYZ> pot_eff_3S1_2x2[2*2];
   
   bool set_args();
   
   void  input_corrdata_from_results();
   void  input_corrdata_from_bin    ();
   void output_corrdata_to_bin      ();
   void output_corrdata_to_text     ();
   
   void  input_NBSdata_from_bin     (const int);
   
   void potential_calculation       (const int);
   void output_XYZdata_to_text      (const int);
   void output_XYZdata_to_bin       (const int);
   
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main() {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   if (input_from_bin_flg) input_corrdata_from_bin();
   else {
      input_corrdata_from_results();
      output_corrdata_to_bin();
   }
   
   if (output_corr_flg) output_corrdata_to_text();
   
   for (int itime=min_time_pot; itime<=max_time_pot; itime++) {
      printf("@@@ Potential calculation at t=%d START\n", itime);
      
      input_NBSdata_from_bin(itime);
      
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
   bool set_args() {
      if      (quark == "c") {
         tmp_mas_3x3_1[0] = mass_nucl; tmp_mas_3x3_1[1] = mass_nucl; tmp_mas_3x3_1[2] = mass_nucl;
         tmp_mas_3x3_2[0] = mass_Lamc; tmp_mas_3x3_2[1] = mass_Sigc; tmp_mas_3x3_2[2] = mass_SigcS;
         
         tmp_mas_2x2_1[0] = mass_nucl; tmp_mas_2x2_1[1] = mass_nucl;
         tmp_mas_2x2_2[0] = mass_Sigc; tmp_mas_2x2_2[1] = mass_SigcS;
         
         ZP_3x3[0] = sq_ZP_nucl*sq_ZP_Lamc; ZP_3x3[1] = sq_ZP_nucl*sq_ZP_Sigc; ZP_3x3[2] = sq_ZP_nucl*sq_ZP_SigcS;
         ZP_2x2[0] = sq_ZP_nucl*sq_ZP_Sigc; ZP_2x2[1] = sq_ZP_nucl*sq_ZP_SigcS;
      }
      else if (quark == "s") {
         tmp_mas_3x3_1[0] = mass_nucl; tmp_mas_3x3_1[1] = mass_nucl; tmp_mas_3x3_1[2] = mass_nucl;
         tmp_mas_3x3_2[0] = mass_Lam ; tmp_mas_3x3_2[1] = mass_Sig ; tmp_mas_3x3_2[2] = mass_SigS;
         
         tmp_mas_2x2_1[0] = mass_nucl; tmp_mas_2x2_1[1] = mass_nucl;
         tmp_mas_2x2_2[0] = mass_Sig ; tmp_mas_2x2_2[1] = mass_SigS;
         
         ZP_3x3[0] = sq_ZP_nucl*sq_ZP_Lam; ZP_3x3[1] = sq_ZP_nucl*sq_ZP_Sig; ZP_3x3[2] = sq_ZP_nucl*sq_ZP_SigS;
         ZP_2x2[0] = sq_ZP_nucl*sq_ZP_Sig; ZP_2x2[1] = sq_ZP_nucl*sq_ZP_SigS;
      }
      else {
         printf("ERROR quark, %s", quark.c_str());
         return false;
      }
      return true;
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_results() {
      conf.set(fconf_list);
      for (int ihad=0; ihad<2; ihad++) {
         for (int ich=0; ich<3; ich++) corr_3x3[ich][ihad].mem_alloc(conf.Nlist());
         for (int ich=0; ich<2; ich++) corr_2x2[ich][ihad].mem_alloc(conf.Nlist());
      }
      printf("Reading Correlators...\n");
      for (int ihad=0; ihad<2; ihad++) for (int i=0; i<conf.Nlist(); i++) {
         for (int ich=0; ich<3; ich++) {
            corr_3x3[ich][ihad](i).mem_alloc(Tsize);
            corr_3x3[ich][ihad](i).input_data_corr(path_corr(ibase_corr, conf(i),
                                                             ihad_3x3_str[ich][ihad], src_shift_str), true);
         }
         for (int ich=0; ich<2; ich++) {
            corr_2x2[ich][ihad](i).mem_alloc(Tsize);
            corr_2x2[ich][ihad](i).input_data_corr(path_corr(ibase_corr, conf(i),
                                                             ihad_2x2_str[ich][ihad], src_shift_str), true);
         }
      }
      printf("Reading Correlators... end\n");
      
      printf("Make Jack-knife samples... \n");
      for (int ihad=0; ihad<2; ihad++) {
         for (int ich=0; ich<3; ich++) corr_3x3[ich][ihad].make_JK_sample(bin_size);
         for (int ich=0; ich<2; ich++) corr_2x2[ich][ihad].make_JK_sample(bin_size);
      }
      printf("Make Jack-knife samples... end\n");
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_bin() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      printf("Reading Correlators...\n");
      for (int ihad=0; ihad<2; ihad++) {
         for (int ich=0; ich<3; ich++)
            corr_3x3[ich][ihad].input_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_3x3_str[ich][ihad]+
                                               ".bin_size"+Bsize_str);
         for (int ich=0; ich<2; ich++)
            corr_2x2[ich][ihad].input_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_2x2_str[ich][ihad]+
                                               ".bin_size"+Bsize_str);
      }
      printf("Reading Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_corrdata_to_bin() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      printf("Writing Correlators...\n");
      for (int ihad=0; ihad<2; ihad++) {
         for (int ich=0; ich<3; ich++)
            corr_3x3[ich][ihad].output_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_3x3_str[ich][ihad]+
                                                ".bin_size"+Bsize_str);
         for (int ich=0; ich<2; ich++)
            corr_2x2[ich][ihad].output_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_2x2_str[ich][ihad]+
                                                ".bin_size"+Bsize_str);
      }
      printf("Writing Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_corrdata_to_text() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      printf("Writing Correlators...\n");
      for (int ihad=0; ihad<2; ihad++) {
         for (int ich=0; ich<3; ich++)
            corr_3x3[ich][ihad].output_data_err(obase+"/"+odir_corr+"/correlator."+ihad_3x3_str[ich][ihad]+
                                                ".bin_size"+Bsize_str+".gnu", 0, true);
         for (int ich=0; ich<2; ich++)
            corr_2x2[ich][ihad].output_data_err(obase+"/"+odir_corr+"/correlator."+ihad_2x2_str[ich][ihad]+
                                                ".bin_size"+Bsize_str+".gnu", 0, true);
      }
      printf("Writing Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void input_NBSdata_from_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Reading NBSwave...\n");
      for (int dt=-1; dt<=+1; dt++) {
         char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
         string time_dt_str(time_dt_c);
         
         for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
            int  ijch = jch + 3*ich;
            wave_3S1_3x3[dt+1][ijch].input_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_3x3_str[ich]+"_"+
                                                    ch_3x3_str[jch]+".spin_1_.t"+time_dt_str+
                                                    ".bin_size"+Bsize_str);
         }
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            int  ijch = jch + 2*ich;
            wave_3S1_2x2[dt+1][ijch].input_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_2x2_str[ich]+"_"+
                                                    ch_2x2_str[jch]+".spin_1_.t"+time_dt_str+
                                                    ".bin_size"+Bsize_str);
         }
      }
      printf("Reading NBSwave... end\n");
   }
   ////// ------------------------------------------------------ //////
   void potential_calculation(const int itime) {
      for (int ijch=0; ijch<3*3; ijch++) pot_eff_3S1_3x3[ijch].mem_alloc(wave_3S1_3x3[0][ijch].Ndata());
      for (int ijch=0; ijch<2*2; ijch++) pot_eff_3S1_2x2[ijch].mem_alloc(wave_3S1_2x2[0][ijch].Ndata());
      
      printf("Calculating potentials... \n");
      
#pragma omp parallel for
      for (int i=0; i<wave_3S1_3x3[0][0].Ndata(); i++){
         ComplexField_XYZ tmp_pot_3x3[3*3], tmp_Rcor_3x3[3][3*3];
         ComplexField_XYZ tmp_pot_2x2[2*2], tmp_Rcor_2x2[3][2*2];
         cdouble          tmp_fac_3x3[3*3], tmp_fac_2x2 [2*2];
         
         for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
            int  ijch = jch + 3*ich;
            tmp_fac_3x3[ijch] = ((ZP_3x3[ich] * corr_3x3[jch][0](i)(itime) * corr_3x3[jch][1](i)(itime)) /
                                 (ZP_3x3[jch] * corr_3x3[ich][0](i)(itime) * corr_3x3[ich][1](i)(itime)));
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor_3x3[dt+1][ijch] = (wave_3S1_3x3[dt+1][ijch](i).rot_proj(ROT_REP_A1) /
                                           (corr_3x3[ich][0](i)(itime+dt) * corr_3x3[ich][1](i)(itime+dt)));
         }
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            int  ijch = jch + 2*ich;
            tmp_fac_2x2[ijch] = ((ZP_2x2[ich] * corr_2x2[jch][0](i)(itime) * corr_2x2[jch][1](i)(itime)) /
                                 (ZP_2x2[jch] * corr_2x2[ich][0](i)(itime) * corr_2x2[ich][1](i)(itime)));
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor_2x2[dt+1][ijch] = (wave_3S1_2x2[dt+1][ijch](i).rot_proj(ROT_REP_A1) /
                                           (corr_2x2[ich][0](i)(itime+dt) * corr_2x2[ich][1](i)(itime+dt)));
         }
         
         Potential::calc_CCpotential_T2(tmp_pot_3x3, tmp_Rcor_3x3[0], tmp_Rcor_3x3[1], tmp_Rcor_3x3[2],
                                        tmp_mas_3x3_1, tmp_mas_3x3_2, tmp_fac_3x3, 3);
         Potential::calc_CCpotential_T2(tmp_pot_2x2, tmp_Rcor_2x2[0], tmp_Rcor_2x2[1], tmp_Rcor_2x2[2],
                                        tmp_mas_2x2_1, tmp_mas_2x2_2, tmp_fac_2x2, 2);
         
         for (int ijch=0; ijch<3*3; ijch++) pot_eff_3S1_3x3[ijch](i) = tmp_pot_3x3[ijch];
         for (int ijch=0; ijch<2*2; ijch++) pot_eff_3S1_2x2[ijch](i) = tmp_pot_2x2[ijch];
         
         printf("Calculating potentials... end: conf=%d\n", i);
      }
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_text(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      printf("Output potentials... \n");
      for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
         int  ijch = jch + 3*ich;
         pot_eff_3S1_3x3[ijch].output_data_err(obase+"/"+odir_pot+"/Pot."+ch_3x3_str[ich]+"_"+
                                               ch_3x3_str[jch]+".Veff3S1.t"+time_str+".bin_size"+
                                               Bsize_str+".gnu", lattice_spacing, true);
      }
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         int  ijch = jch + 2*ich;
         pot_eff_3S1_2x2[ijch].output_data_err(obase+"/"+odir_pot+"/Pot."+ch_2x2_str[ich]+"_"+
                                               ch_2x2_str[jch]+".Veff3S1.t"+time_str+".bin_size"+
                                               Bsize_str+".gnu", lattice_spacing, true);
      }
      printf("Output potentials... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      printf("Output potentials... \n");
      for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
         int  ijch = jch + 3*ich;
         pot_eff_3S1_3x3[ijch].output_data_bin_reduce(obase+"/"+odir_pot+"/Pot."+ch_3x3_str[ich]+"_"+
                                                      ch_3x3_str[jch]+".Veff3S1.t"+time_str+".bin_size"+
                                                      Bsize_str+".bin", lattice_spacing, false);
      }
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         int  ijch = jch + 2*ich;
         pot_eff_3S1_2x2[ijch].output_data_bin_reduce(obase+"/"+odir_pot+"/Pot."+ch_2x2_str[ich]+"_"+
                                                      ch_2x2_str[jch]+".Veff3S1.t"+time_str+".bin_size"+
                                                      Bsize_str+".bin", lattice_spacing, false);
      }
      printf("Output potentials... end\n");
   }
}
