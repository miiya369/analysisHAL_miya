#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int    Tsize = 64;
   int    Lsize = 32;
   
   int bin_size = 57;
   
   int min_time_wave = 9;
   int max_time_wave = 13;
   
   string quark = "c";
   
   string wdir       = "/gwfefs/data/G18002/miyamoto/work.pacs-cs/bridge-1.4.3-hal.v01/src2/run_rela";
   string ibase      = wdir +"/results_"+ quark;
   string ibase_corr = wdir +"/results.split/"+ quark +"/results.split.ave";
   string obase      = wdir +"/analysis/"+ quark +"/bin";
   string fconf_list = wdir +"/conf.lst";
   
   string path_wave(string dir_base, string dir_channel, string conf, int time, string shift, string footer) {
      char tmp_c[8]; snprintf(tmp_c, sizeof(tmp_c), "%+04d", time);
      return (dir_base +"/"+ dir_channel +"/"+ conf +"/NBSwave."+ tmp_c +"+"+ shift +"."+ conf +"."+ footer);
   }
   string path_corr(string dir_base, string conf, string hadron, string shift) {
      return (dir_base +"/correlator.PS.dir/"+ conf + "/" + hadron +"_correlator.+"+ shift +"."+ conf);
   }
   
   string ichdir_str[5] = {"Prot__Lamb_", "Prot__SigZ_", "Neut__SigP_", "SigSZ_Prot_", "SigSP_Neut_"};
   int    identif_ch[5] = {HH_OctOct    , HH_OctOct    , HH_OctOct    , HH_DecOct    , HH_DecOct    };
   
   string ihad_str[4] = {"proton_CG05_CG05", "Lambda_CG05_CG05", "Sigma_CG05_CG05", "Sigma32"};
   int    Nhad        = 4;
   
   string ch3_3x3_str[3] = {"L"+quark+"N12_", "S"+quark+"N12_", "S"+quark+"SN12"};
   string ch3_2x2_str[2] = {"S"+quark+"N32_", "S"+quark+"SN32"};
   
   string src_shift_str  = "A64.000.000.000";
   
   cmatrix Pmat_5x5(5);
   cmatrix Pmat_3x3(3);
   
   bool set_args();
   
   void io_corr    (const string);
   void io_NBS_all (const int, const int);
   void io_NBS_spin(const int, const int);
   
   void input_NBS_sub(ComplexField_AXYZB&, const int, const int, const int, const string);
   
   void spprj_NBS_sub(const ComplexField_AXYZB&, ComplexFieldMatrix<ComplexField_XYZ> wave_sprj[2],
                      const int, const int);
   template <class X>
   void isprj_NBS_sub(STATISTICS<X>*, STATISTICS<X>*, ComplexFieldMatrix<X>&, const int, const int);
   
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main() {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   for (int ihad=0; ihad<Nhad; ihad++) io_corr(ihad_str[ihad]);
   
   for (int itime=max_time_wave; itime>=min_time_wave; itime--)
   {
      //io_NBS_all(itime, 2);
      io_NBS_all(itime, 3);
      
      //io_NBS_spin(itime, 2);
      io_NBS_spin(itime, 3);
   }
   
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args() {
      cmatrix Pmat(5);
      // Isospin projection matrix for S=-1 channel
      Pmat(0,0) = 1.0; Pmat(0,1) =           0.0 ; Pmat(0,2) =          0.0 ; Pmat(0,3) =           0.0 ; Pmat(0,4) =          0.0 ;
      Pmat(1,0) = 0.0; Pmat(1,1) = -sqrt(1.0/3.0); Pmat(1,2) = sqrt(2.0/3.0); Pmat(1,3) =           0.0 ; Pmat(1,4) =          0.0 ;
      Pmat(2,0) = 0.0; Pmat(2,1) =  sqrt(2.0/3.0); Pmat(2,2) = sqrt(1.0/3.0); Pmat(2,3) =           0.0 ; Pmat(2,4) =          0.0 ;
      Pmat(3,0) = 0.0; Pmat(3,1) =           0.0 ; Pmat(3,2) =          0.0 ; Pmat(3,3) = -sqrt(1.0/3.0); Pmat(3,4) = sqrt(2.0/3.0);
      Pmat(4,0) = 0.0; Pmat(4,1) =           0.0 ; Pmat(4,2) =          0.0 ; Pmat(4,3) =  sqrt(2.0/3.0); Pmat(4,4) = sqrt(1.0/3.0);
      
      for (int i=0; i<5; i++) for (int j=0; j<5; j++) Pmat_5x5(i,j) = Pmat(i,j);
      for (int i=0; i<3; i++) for (int j=0; j<3; j++) Pmat_3x3(i,j) = Pmat(i,j);
      
      return true;
   }
   ////// ------------------------------------------------------ //////
   void io_corr(const string had_name) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      
      anaHAL::NameList conf(fconf_list);
      STATISTICS<ComplexField_T> corr(conf.Nlist());
      
      printf("Reading Correlator[%s]... ", had_name.c_str()); fflush(stdout);
      for (int i=0; i<conf.Nlist(); i++) {
         corr(i).mem_alloc(Tsize);
         corr(i).input_data_corr(path_corr(ibase_corr, conf(i), had_name, src_shift_str), true);
      }
      printf("end\n");
      
      printf("Make Jack-knife samples... "); fflush(stdout);
      corr.make_JK_sample(bin_size);
      printf("end\n");
      
      printf("Writing Correlator[%s]... ", had_name.c_str()); fflush(stdout);
      corr.output_data_bin(obase+"/correlator."+had_name+".bin_size"+Bsize_str);
      printf("end\n");
   }
   ////// ------------------------------------------------------ //////
   void io_NBS_all(const int itime, const int Nch) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      if (Nch < 2 || 3 < Nch) ERROR_COMMENTS("Invalid #.channel.");
      
      anaHAL::NameList conf(fconf_list);
      STATISTICS<ComplexField_AXYZB> *wave1 = new STATISTICS<ComplexField_AXYZB> [ Nch   * Nch   ];
      STATISTICS<ComplexField_AXYZB> *wave2 = new STATISTICS<ComplexField_AXYZB> [(Nch-1)*(Nch-1)];
      
      for (int ijch=0; ijch< Nch   * Nch   ; ijch++) wave1[ijch].mem_alloc(conf.Nlist());
      for (int ijch=0; ijch<(Nch-1)*(Nch-1); ijch++) wave2[ijch].mem_alloc(conf.Nlist());
      
      ComplexField_AXYZB tmp_wave;
      ComplexFieldMatrix<ComplexField_AXYZB> wave_mat(2*Nch-1);
      
      printf("Reading NBSwave[it=%d]...\n", itime);
      for (int i=0; i<conf.Nlist(); i++) {
         for (int ich=0; ich<2*Nch-1; ich++) for (int jch=0; jch<2*Nch-1; jch++) {
            input_NBS_sub(tmp_wave, ich, jch, itime, conf(i));
            wave_mat(ich,jch) = tmp_wave;
         }
         isprj_NBS_sub(wave1, wave2, wave_mat, Nch, i);
         
         printf("Reading NBSwave[it=%d]... end: conf=%d\n", itime, i);
      }
      printf("Make Jack-knife samples... "); fflush(stdout);
      for (int ijch=0; ijch< Nch   * Nch   ; ijch++) wave1[ijch].make_JK_sample(bin_size);
      for (int ijch=0; ijch<(Nch-1)*(Nch-1); ijch++) wave2[ijch].make_JK_sample(bin_size);
      printf("end\n");
      
      printf("Writing NBSwave[it=%d]... ", itime); fflush(stdout);
      for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
         int ijch = jch + Nch*ich;
         wave1[ijch].output_data_bin(obase+"/NBSwave."+ch3_3x3_str[ich]+"_"+ch3_3x3_str[jch]+
                                     ".allspin.t"+time_str+".bin_size"+Bsize_str);
      }
      for (int ich=0; ich<Nch-1; ich++) for (int jch=0; jch<Nch-1; jch++) {
         int ijch = jch + (Nch-1)*ich;
         wave2[ijch].output_data_bin(obase+"/NBSwave."+ch3_2x2_str[ich]+"_"+ch3_2x2_str[jch]+
                                     ".allspin.t"+time_str+".bin_size"+Bsize_str);
      }
      printf("end\n");
      
      delete [] wave1;
      delete [] wave2;
   }
   ////// ------------------------------------------------------ //////
   void io_NBS_spin(const int itime, const int Nch) {
      char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",    itime); string  time_str( time_c);
      
      if (Nch < 2 || 3 < Nch) ERROR_COMMENTS("Invalid #.channel.");
      
      anaHAL::NameList conf(fconf_list);
      STATISTICS<ComplexField_XYZ> *wave1[2], *wave2[2]; // [2]=spin 0 (2 for DecOct-DecOct) and 1
      
      for (int ispin=0; ispin<2; ispin++) {
         wave1[ispin] = new STATISTICS<ComplexField_XYZ> [ Nch   * Nch   ];
         wave2[ispin] = new STATISTICS<ComplexField_XYZ> [(Nch-1)*(Nch-1)];
         
         for (int ijch=0; ijch< Nch   * Nch   ; ijch++) wave1[ispin][ijch].mem_alloc(conf.Nlist());
         for (int ijch=0; ijch<(Nch-1)*(Nch-1); ijch++) wave2[ispin][ijch].mem_alloc(conf.Nlist());
      }
      
      ComplexField_AXYZB tmp_wave;
      ComplexFieldMatrix<ComplexField_XYZ> wave_mat[2]; // [2]=spin, 0 and 1
      for (int ispin=0; ispin<2; ispin++)  wave_mat[ispin].init(2*Nch-1);
      
      printf("Reading NBSwave[it=%d]...\n", itime);
      for (int i=0; i<conf.Nlist(); i++) {
         for (int ich=0; ich<2*Nch-1; ich++) for (int jch=0; jch<2*Nch-1; jch++)
         {
            input_NBS_sub(tmp_wave, ich, jch, itime, conf(i));
            spprj_NBS_sub(tmp_wave, wave_mat, ich, jch);
         }
         for (int ispin=0; ispin<2; ispin++)
            isprj_NBS_sub(wave1[ispin], wave2[ispin], wave_mat[ispin], Nch, i);
         
         printf("Reading NBSwave[it=%d]... end: conf=%d\n", itime, i);
      }
      printf("Make Jack-knife samples... "); fflush(stdout);
      for (int ispin=0; ispin<2; ispin++) {
         for (int ijch=0; ijch< Nch   * Nch   ; ijch++) wave1[ispin][ijch].make_JK_sample(bin_size);
         for (int ijch=0; ijch<(Nch-1)*(Nch-1); ijch++) wave2[ispin][ijch].make_JK_sample(bin_size);
      }
      printf("end\n");
      
      printf("Writing NBSwave[it=%d]... ", itime); fflush(stdout);
      for (int ispin=0; ispin<2; ispin++) {
         char spin_c[16];
         snprintf(spin_c, sizeof(spin_c), ".spin_%d_.t", ispin); string spin_str(spin_c);
         
         for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
            int ijch = jch + Nch*ich;
            
            if (ispin == 0 && identif_ch[ich] == HH_DecOct)
               snprintf(spin_c, sizeof(spin_c), ".spin_2_.t"); spin_str = spin_c;
            
            wave1[ispin][ijch].output_data_bin(obase+"/NBSwave."+ch3_3x3_str[ich]+"_"+ch3_3x3_str[jch]+
                                               spin_str+time_str+".bin_size"+Bsize_str);
         }
         for (int ich=0; ich<Nch-1; ich++) for (int jch=0; jch<Nch-1; jch++) {
            int ijch = jch + (Nch-1)*ich;
            
            if (ispin == 0 && identif_ch[ich] == HH_DecOct)
               snprintf(spin_c, sizeof(spin_c), ".spin_2_.t"); spin_str = spin_c;
            
            wave2[ispin][ijch].output_data_bin(obase+"/NBSwave."+ch3_2x2_str[ich]+"_"+ch3_2x2_str[jch]+
                                               spin_str+time_str+".bin_size"+Bsize_str);
         }
      }
      printf("end\n");
      
      for (int ispin=0; ispin<2; ispin++) {
         delete [] wave1[ispin];
         delete [] wave2[ispin];
      }
   }
   ////// ------------------------------------------------------ //////
   void input_NBS_sub(ComplexField_AXYZB& wave, const int ich, const int jch, const int itime, const string conf)
   {
      ComplexField_AXYZB wave_in;
      string             ch_str, ch_footer;
      int                factor = 1;
      
      if      (identif_ch[ich] == HH_OctOct && identif_ch[jch] == HH_OctOct) {
         ch_str    = "2Boct";
         ch_footer = "oo_CG05.oo_CG05";
         factor    = +1;
         wave_in.mem_alloc(2*2  , Lsize, 2*2  );
      }
      else if (identif_ch[ich] == HH_OctOct && identif_ch[jch] == HH_DecOct) {
         ch_str    = "ooxdo";
         ch_footer = "oo_CG05.do_CG0K";
         factor    = +1;
         wave_in.mem_alloc(2*2  , Lsize, 2*2*3);
      }
      else if (identif_ch[ich] == HH_DecOct && identif_ch[jch] == HH_OctOct) {
         ch_str    = "doxoo";
         ch_footer = "do_CG0K.oo_CG05";
         factor    = -1;
         wave_in.mem_alloc(2*2*3, Lsize, 2*2  );
      }
      else if (identif_ch[ich] == HH_DecOct && identif_ch[jch] == HH_DecOct) {
         ch_str    = "2Bdow";
         ch_footer = "do_CG0K.do_CG0K";
         factor    = -1;
         wave_in.mem_alloc(2*2*3, Lsize, 2*2*3);
      }
      else ERROR_COMMENTS("Invalid channel.");
      
      wave_in.input_data_bin(path_wave(ibase, "NBS_"+ch_str+"."+ichdir_str[ich]+"_"+
                                       ichdir_str[jch]+".dir", conf, +itime, "Ave.000.000.000", ch_footer));
      wave  = wave_in;
      wave_in.input_data_bin(path_wave(ibase, "NBS_"+ch_str+"."+ichdir_str[ich]+"_"+
                                       ichdir_str[jch]+".dir", conf, -itime, "Ave.000.000.000", ch_footer));
      wave += wave_in * factor;
      wave_in.input_data_bin(path_wave(ibase, "NBS_"+ch_str+"."+ichdir_str[ich]+"_"+
                                       ichdir_str[jch]+".dir", conf, +itime, "Av2.000.000.000", ch_footer));
      wave += wave_in;
      wave_in.input_data_bin(path_wave(ibase, "NBS_"+ch_str+"."+ichdir_str[ich]+"_"+
                                       ichdir_str[jch]+".dir", conf, -itime, "Av2.000.000.000", ch_footer));
      wave += wave_in * factor;
      wave /= 4.0;
   }
   ////// ------------------------------------------------------ //////
   void spprj_NBS_sub(const ComplexField_AXYZB& wave_org, ComplexFieldMatrix<ComplexField_XYZ> wave_sprj[2],
                      const int ich, const int jch)
   {
      if      (identif_ch[ich] == HH_OctOct && identif_ch[jch] == HH_OctOct)
         wave_sprj[0](ich,jch)  = wave_org.spin_proj(identif_ch[ich], 0, 0, identif_ch[jch], 0, 0);
      
      else if (identif_ch[ich] == HH_OctOct && identif_ch[jch] == HH_DecOct)
         wave_sprj[0](ich,jch)  = wave_org.spin_proj(identif_ch[ich], 0, 0, identif_ch[jch], 2, 0);
      
      else if (identif_ch[ich] == HH_DecOct && identif_ch[jch] == HH_OctOct)
         wave_sprj[0](ich,jch)  = wave_org.spin_proj(identif_ch[ich], 2, 0, identif_ch[jch], 0, 0);
      
      else if (identif_ch[ich] == HH_DecOct && identif_ch[jch] == HH_DecOct) {
         wave_sprj[0](ich,jch)  = wave_org.spin_proj(identif_ch[ich], 2,+2, identif_ch[jch], 2,+2);
         wave_sprj[0](ich,jch) += wave_org.spin_proj(identif_ch[ich], 2,+1, identif_ch[jch], 2,+1);
         wave_sprj[0](ich,jch) += wave_org.spin_proj(identif_ch[ich], 2, 0, identif_ch[jch], 2, 0);
         wave_sprj[0](ich,jch) += wave_org.spin_proj(identif_ch[ich], 2,-1, identif_ch[jch], 2,-1);
         wave_sprj[0](ich,jch) += wave_org.spin_proj(identif_ch[ich], 2,-2, identif_ch[jch], 2,-2);
         wave_sprj[0](ich,jch) /= 5.0;
      }
      
      wave_sprj[1](ich,jch)  = wave_org.spin_proj(identif_ch[ich], 1,+1, identif_ch[jch], 1,+1);
      wave_sprj[1](ich,jch) += wave_org.spin_proj(identif_ch[ich], 1, 0, identif_ch[jch], 1, 0);
      wave_sprj[1](ich,jch) += wave_org.spin_proj(identif_ch[ich], 1,-1, identif_ch[jch], 1,-1);
      wave_sprj[1](ich,jch) /= 3.0;
   }
   ////// ------------------------------------------------------ //////
   template <class X>
   void isprj_NBS_sub(STATISTICS<X> *w1, STATISTICS<X> *w2, ComplexFieldMatrix<X>& wave_mat,
                      const int Nch, const int i)
   {
      if (Nch == 2) {
         wave_mat = Pmat_3x3 * wave_mat * Pmat_3x3.T();
         
         w1[0](i) = wave_mat(0,0);
         w1[1](i) = wave_mat(0,1);
         w1[2](i) = wave_mat(1,0);
         w1[3](i) = wave_mat(1,1);
         w2[0](i) = wave_mat(2,2);
      }
      else if (Nch == 3) {
         wave_mat = Pmat_5x5 * wave_mat * Pmat_5x5.T();
         
         w1[0](i) = wave_mat(0,0);
         w1[1](i) = wave_mat(0,1);
         w1[2](i) = wave_mat(0,3);
         w1[3](i) = wave_mat(1,0);
         w1[4](i) = wave_mat(1,1);
         w1[5](i) = wave_mat(1,3);
         w1[6](i) = wave_mat(3,0);
         w1[7](i) = wave_mat(3,1);
         w1[8](i) = wave_mat(3,3);
         
         w2[0](i) = wave_mat(2,2);
         w2[1](i) = wave_mat(2,4);
         w2[2](i) = wave_mat(4,2);
         w2[3](i) = wave_mat(4,4);
      }
      else ERROR_COMMENTS("Invalid #.channel.");
   }
}
