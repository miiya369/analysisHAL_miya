#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>
#include <Potential.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int    Tsize = 32;
   int    Lsize = 16;
   double lattice_spacing = 0.1215;
   
   int bin_size = 1;
   
   int min_time_pot = 3;
   int max_time_pot = 12;
   
   bool pre_read_flg    = true;
   bool output_corr_flg = true;
   bool output_pot_flg  = true;
   
   string alpha[2]   = {"000", "016"};
   
   string dir_base   = "/nas/nasbee01/miyamoto/Data/su3/alpha_src";
   string ibase_corr = dir_base + "/results.Alpha.";
   string ibase_nbs  = dir_base + "/results.Alpha.";
   string obase      = dir_base + "/analysis";
   string fconf_list = dir_base + "/conf.lst";
   
   string odir_bin   = "bin";
   string odir_corr  = "corr";
   string odir_pot   = "pot";
   string odir_wave  = "wave";
   
   string path_wave(string dir_base, string dir_channel, string conf, int time, string shift, string footer) {
      char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "%+04d", time);
      return (dir_base +"/"+ dir_channel +"/"+ conf +"/NBSwave."+ tmp_c +"+"+ shift +"."+ conf +"."+ footer);
   }
   string path_corr(string dir_base, string conf, string hadron, string shift) {
      return (dir_base +"/correlator.PS.dir/"+ conf + "/" + hadron +"_correlator.+"+ shift +"."+ conf);
   }
   
   string ihad_str  [2] = {"proton_NR_NR", "OMEGA"};
   
   string ch_nbs_str[2] = {"NN0_DD0", "DD0_DD0"};
   string ch_pot_str[2] = {"NN0", "DD0"};
   
   string src_shift_str = "000.000.000.000";
   
   double mass_nuc_a00 = 1.3;
   double mass_nuc_a08 = 1.3;
   double mass_nuc_a16 = 1.3;
   
   double mass_del_a00 = 1.4;
   double mass_del_a08 = 1.4;
   double mass_del_a16 = 1.4;
   
   double sqZP_nuc_a00 = 0.27;
   double sqZP_nuc_a08 = 0.27;
   double sqZP_nuc_a16 = 0.27;
   
   double sqZP_del_a00 = 0.012;
   double sqZP_del_a08 = 0.012;
   double sqZP_del_a16 = 0.012;
   
   double  had_mas[2][2]; // [channel][alpha]
   cdouble Zfactor[2][2];
   
   anaHAL::NameList conf;
   STATISTICS<ComplexField_T  > corr[2][2];  // [channel][alpha]
   STATISTICS<ComplexField_XYZ> wave[3][2*2];
   STATISTICS<ComplexField_XYZ>  pot[2*2], detR;
   
   bool set_args();
   
   void  input_corrdata_from_results();
   void  input_corrdata_from_bin    ();
   void output_corrdata_to_bin      ();
   void output_corrdata_to_text     ();
   
   void io_NBSdata_from_results(const int);
   void input_NBSdata_from_bin (const int);
   
   void potential_calculation  (const int);
   void output_XYZdata_to_text (const int);
   void output_XYZdata_to_bin  (const int);
   
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main(int argc, char **argv) {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   if (pre_read_flg) {
      input_corrdata_from_results();
      output_corrdata_to_bin();
      
      for (int itime=min_time_pot-1; itime<=max_time_pot+1; itime++)
         io_NBSdata_from_results(itime);
   }
   else {
      input_corrdata_from_bin();
      if (output_corr_flg) output_corrdata_to_text();
      
      for (int itime=min_time_pot; itime<=max_time_pot; itime++) {
         printf("@@@ Potential calculation at t=%d START\n", itime);
         
         input_NBSdata_from_bin(itime);
         
         if (output_pot_flg) {
            potential_calculation (itime);
            output_XYZdata_to_text(itime);
            output_XYZdata_to_bin (itime);
         }
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args() {
      had_mas[0][0] = mass_nuc_a00;
      had_mas[0][1] = mass_nuc_a16;
      
      had_mas[1][0] = mass_del_a00;
      had_mas[1][1] = mass_del_a16;
      
      Zfactor[0][0] = pow(sqZP_nuc_a00, 2);
      Zfactor[0][1] = pow(sqZP_nuc_a16, 2);
      Zfactor[1][0] = pow(sqZP_del_a00, 2);
      Zfactor[1][1] = pow(sqZP_del_a16, 2);
      
      return true;
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_results() {
      conf.set(fconf_list);
      
      printf("Reading Correlators...\n");
      for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++) {
         corr[ich][ia].mem_alloc(conf.Nlist());
         for (int i=0; i<conf.Nlist(); i++) {
            corr[ich][ia](i).mem_alloc(Tsize);
            corr[ich][ia](i).input_data_corr(path_corr(ibase_corr+alpha[ia], conf(i),
                                                       ihad_str[ich], src_shift_str), true);
         }
      }
      printf("Reading Correlators... end\n");
      
      printf("Make Jack-knife samples... \n");
      for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++)
         corr[ich][ia].make_JK_sample(bin_size);
      printf("Make Jack-knife samples... end\n");
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_bin() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Reading Correlators...\n");
      for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++)
         corr[ich][ia].input_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_str[ich]+
                                      ".alpha"+alpha[ia]+".bin_size"+Bsize_str);
      printf("Reading Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_corrdata_to_bin() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Writing Correlators...\n");
      for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++)
         corr[ich][ia].output_data_bin(obase+"/"+odir_bin+"/correlator."+ihad_str[ich]+
                                       ".alpha"+alpha[ia]+".bin_size"+Bsize_str);
      printf("Writing Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_corrdata_to_text() {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Writing Correlators...\n");
      for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++)
         corr[ich][ia].output_data_err(obase+"/"+odir_corr+"/correlator."+ihad_str[ich]+
                                       ".alpha"+alpha[ia]+".bin_size"+Bsize_str+".gnu", 0, true);
      printf("Writing Correlators... end\n");
   }
   ////// ------------------------------------------------------ //////
   void io_NBSdata_from_results(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      conf.set(fconf_list);
      
      for (int iach=0; iach<2*2; iach++) wave[0][iach].mem_alloc(conf.Nlist());
      
      ComplexField_XYZ tmp_wave[4]; for (int i=0; i<4; i++) tmp_wave[i].mem_alloc(Lsize);
      ComplexField_XYZ  in_wave(Lsize);
      string ich_nbs_str[2] = {"DelPP_DelM_", "DelP__DelZ_"};
      
      printf("Reading & Writing NBSwave[it=%d]...\n", itime);
      for (int ia=0; ia<2; ia++) {
         for (int i=0; i<conf.Nlist(); i++) {
	    for (int ich=0; ich<2; ich++) {
	       in_wave.input_data_bin(path_wave(ibase_nbs+alpha[ia],
						"Proj.NBS_ooxdd.Prot__Neut__"+ich_nbs_str[ich]+
						".dir/spin1z+0.3z+0",
						conf(i),  itime, src_shift_str, "oo_NR.dd_NR"));
	       tmp_wave[ich]  = in_wave;
	       
	       in_wave.input_data_bin(path_wave(ibase_nbs+alpha[ia],
						"Proj.NBS_ooxdd.Prot__Neut__"+ich_nbs_str[ich]+
						".dir/spin1z+0.3z+0",
						conf(i), -itime, src_shift_str, "oo_NR.dd_NR"));
	       tmp_wave[ich] += in_wave;
	       tmp_wave[ich] /= 2.0;
	    }
            wave[0][ia+2*0](i) = (tmp_wave[0] - tmp_wave[1]) / 2.0;
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int ijch = jch + 2*ich;
	       in_wave.input_data_bin(path_wave(ibase_nbs+alpha[ia],
						"Proj.NBS_2Bdec."+ich_nbs_str[ich]+"_"+
						ich_nbs_str[jch]+".dir/spin3z+0.3z+0",
						conf(i),  itime, src_shift_str, "dd_NR.dd_NR"));
	       tmp_wave[ijch]  = in_wave;
	       
               in_wave.input_data_bin(path_wave(ibase_nbs+alpha[ia],
						"Proj.NBS_2Bdec."+ich_nbs_str[ich]+"_"+
						ich_nbs_str[jch]+".dir/spin3z+0.3z+0",
						conf(i), -itime, src_shift_str, "dd_NR.dd_NR"));
	       tmp_wave[ijch] += in_wave;
	       tmp_wave[ijch] /= 2.0;
            }
            wave[0][ia+2*1](i) = (tmp_wave[0] - tmp_wave[1] - tmp_wave[2] + tmp_wave[3]) / 4.0;
         }
         for (int ich=0; ich<2; ich++) wave[0][ia+2*ich].make_JK_sample(bin_size);
         for (int ich=0; ich<2; ich++)
            wave[0][ia+2*ich].output_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_nbs_str[ich]+
                                              ".alpha"+alpha[ia]+".t"+time_str+".bin_size"+Bsize_str);
      } // ia
      printf("Reading & Writing NBSwave[it=%d]... end\n", itime);
   }
   ////// ------------------------------------------------------ //////
   void input_NBSdata_from_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      printf("Reading NBSwave...\n");
      for (int dt=-1; dt<=+1; dt++) {
         char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
         string time_dt_str(time_dt_c);
         
         for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++) {
            int iach = ia + 2*ich;
            wave[dt+1][iach].input_data_bin(obase+"/"+odir_bin+"/NBSwave."+ch_nbs_str[ich]+
                                            ".alpha"+alpha[ia]+".t"+time_dt_str+
                                            ".bin_size"+Bsize_str);
         }
      }
      printf("Reading NBSwave... end\n");
   }
   ////// ------------------------------------------------------ //////
   void potential_calculation(const int itime) {
      for (int ijch=0; ijch<2*2; ijch++) pot[ijch].mem_alloc(wave[0][ijch].Ndata());
      detR.mem_alloc(wave[0][0].Ndata());
      
      printf("Calculating potentials... \n");
      
      #pragma omp parallel for
      for (int i=0; i<wave[0][0].Ndata(); i++) {
         ComplexField_XYZ Rcor[3][2*2], KR[2*2];
         cdouble          fac[2];
         cmatrix          Rmat(2), Kmat(2), Vmat(2);
         
         for (int ia=0; ia<2; ia++)
            fac[ia] = ((Zfactor[0][ia] * pow(corr[1][ia](i)(itime), 2)) /
                       (Zfactor[1][ia] * pow(corr[0][ia](i)(itime), 2)));
         
         for (int dt=-1; dt<=+1; dt++) for (int ia=0; ia<2; ia++) {
            Rcor[dt+1][ia+2*0] = ((wave[dt+1][ia+2*0](i) -
                                   wave[dt+1][ia+2*0](i).rot_proj(ROT_REP_A1)) / sfunc::cfield_Ylm(2,0,Lsize) /
                                  (corr[0][ia](i)(itime+dt)*corr[0][ia](i)(itime+dt))).rot_proj(ROT_REP_A1);
            
            Rcor[dt+1][ia+2*1] = ( wave[dt+1][ia+2*1](i).rot_proj(ROT_REP_A1) * fac[ia] /
                                  (corr[1][ia](i)(itime+dt)*corr[1][ia](i)(itime+dt)));
         }
         for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++) {
            int iach = ia + 2*ich;
            KR[iach] = Potential::get_potential_T2_nume(Rcor[0][iach],
                                                        Rcor[1][iach],
                                                        Rcor[2][iach], had_mas[ich][ia], had_mas[ich][ia]);
            pot[iach](i).mem_alloc(Lsize);
         }
         detR(i).mem_alloc(Lsize);
         
         for (int n=0; n<pot[0](i).data_size(); n++) {
            for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++) {
               int iach = ia + 2*ich;
               Rmat(ich,ia) = Rcor[1][iach](n);
               Kmat(ich,ia) = KR     [iach](n);
            }
            Vmat = Kmat * Rmat.inverce();
            
            for (int ich=0; ich<2; ich++) for (int ia=0; ia<2; ia++) {
               int iach = ia + 2*ich;
               pot[iach](i)(n) = Vmat(ich,ia);
            }
            detR(i)(n) = Rmat.det();
         }
      }
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_text(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      printf("Output potentials... \n");
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         int ijch = jch + 2*ich;
         pot[ijch].output_data_err(obase+"/"+odir_pot+"/Pot."+ch_pot_str[ich]+"_"+ch_pot_str[jch]+
                                   ".a"+alpha[0]+"_"+alpha[1]+".t"+time_str+".bin_size"+Bsize_str+".gnu",
                                   lattice_spacing, true);
      }
      detR.output_data_err(obase+"/"+odir_pot+"/detR.a"+alpha[0]+"_"+alpha[1]+".t"+time_str+".bin_size"+
                           Bsize_str+".gnu", lattice_spacing, true);
      
      for (int ia=0; ia<2; ia++) {
         for (int i=0; i<wave[1][0].Ndata(); i++)
            detR(i) = ((wave[1][ia+2*0](i) -
                        wave[1][ia+2*0](i).rot_proj(ROT_REP_A1)) / sfunc::cfield_Ylm(2,0,Lsize) /
                       (corr[0][ia](i)(itime)*corr[0][ia](i)(itime))).rot_proj(ROT_REP_A1);
         detR.output_data_err(obase+"/"+odir_wave+"/Rcorr."+ch_nbs_str[0]+".alpha"+alpha[ia]+
                              ".t"+time_str+".bin_size"+Bsize_str+".gnu", lattice_spacing, true);
         
         for (int i=0; i<wave[1][0].Ndata(); i++) {
            cdouble fac = ((Zfactor[0][ia] * pow(corr[1][ia](i)(itime), 2)) /
                           (Zfactor[1][ia] * pow(corr[0][ia](i)(itime), 2)));
            detR(i) = ( wave[1][ia+2*1](i).rot_proj(ROT_REP_A1) * fac /
                       (corr[1][ia](i)(itime)*corr[1][ia](i)(itime)));
         }
         detR.output_data_err(obase+"/"+odir_wave+"/Rcorr."+ch_nbs_str[1]+".alpha"+alpha[ia]+
                              ".t"+time_str+".bin_size"+Bsize_str+".gnu", lattice_spacing, true);
      }
      printf("Output potentials... end\n");
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_bin(const int itime) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      printf("Output potentials... \n");
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         int ijch = jch + 2*ich;
         pot[ijch].output_data_bin_reduce(obase+"/"+odir_pot+"/Pot."+ch_pot_str[ich]+"_"+ch_pot_str[jch]+
                                          ".a"+alpha[0]+"_"+alpha[1]+".t"+time_str+".bin_size"+
                                          Bsize_str+".bin", lattice_spacing, false);
      }
      printf("Output potentials... end\n");
   }
}
