#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int    Tsize = 64;
   int    Lsize = 32;
   
   int min_time = 6;
   int max_time = 16;
   
   string ibase = "/home/miiya369/data/pacs-cs";
   string obase = "/home/miiya369/data/pacs-cs/analysis_20180129";
   
   string src_shift_str = "A64.000.000.000";
   
   string path_conf_list(string sc, int ens, int bin_size) {
      char tmp_c[256];
      if      (sc == "s") snprintf(tmp_c, sizeof(tmp_c), "results_kekb_pacs-cs_bridge_Clover/Ens%d_bin%d", ens, bin_size);
      else if (sc == "c") snprintf(tmp_c, sizeof(tmp_c), "results_kekb_pacs-cs_bridge_RHQ/Ens%d_bin%d"   , ens, bin_size);
      else                snprintf(tmp_c, sizeof(tmp_c), "ERROR:_Unknown_sc");
      return (ibase+"/"+tmp_c+"/conf.lst");
   }
   string path_corr_res(string sc, int ens, int bin_size, string conf, string hadron) {
      char tmp_c[256];
      if      (sc == "s") snprintf(tmp_c, sizeof(tmp_c), "results_kekb_pacs-cs_bridge_Clover/Ens%d_bin%d", ens, bin_size);
      else if (sc == "c") snprintf(tmp_c, sizeof(tmp_c), "results_kekb_pacs-cs_bridge_RHQ/Ens%d_bin%d"   , ens, bin_size);
      else                snprintf(tmp_c, sizeof(tmp_c), "ERROR:_Unknown_sc");
      return (ibase+"/"+tmp_c+"/correlator.PS.dir/"+conf+"/"+hadron+"_correlator.+"+src_shift_str+"."+conf);
   }
   string path_wave_res(string sc, int ens, int bin_size, string dir_channel, string conf, int time) {
      char tmp_c[256], tmp_c_t[8]; snprintf(tmp_c_t, sizeof(tmp_c_t), "%+04d", time);
      if      (sc == "s") snprintf(tmp_c, sizeof(tmp_c), "results_kekb_pacs-cs_bridge_Clover/Ens%d_bin%d", ens, bin_size);
      else if (sc == "c") snprintf(tmp_c, sizeof(tmp_c), "results_kekb_pacs-cs_bridge_RHQ/Ens%d_bin%d"   , ens, bin_size);
      else                snprintf(tmp_c, sizeof(tmp_c), "ERROR:_Unknown_sc");
      return (ibase+"/"+tmp_c+"/"+dir_channel+"/"+conf+"/NBSwave."+
              tmp_c_t+"+"+src_shift_str+"."+conf+".NUC_CG05.NUC_CG05");
   }
   string path_corr_bin(string sc, int ens, int bin_size, string hadron) {
      char tmp_c[256];
      snprintf(tmp_c, sizeof(tmp_c), "%s_quark/ens%d/bin/correlator.%s.bin_size%02d",
               sc.c_str(), ens, hadron.c_str(), bin_size);
      return (obase+"/"+tmp_c);
   }
   string path_wave_bin(string sc, int ens, int bin_size, string ch_snk, string ch_src, string spin, int time) {
      char tmp_c[256];
      snprintf(tmp_c, sizeof(tmp_c), "%s_quark/ens%d/bin/NBSwave.%s_%s.%s.t%02d.bin_size%2d",
               sc.c_str(), ens, ch_snk.c_str(), ch_src.c_str(), spin.c_str(), time, bin_size);
      return (obase+"/"+tmp_c);
   }
   
   cmatrix Pmat_LxN (3);
   cmatrix Pmat_LxLx(3);
   
   bool set_args();
   
   void in_out_CORRdata        (const string sc, const int ens, const int bin_size, const string hadron);
   
   void in_out_NBSdata_2x2_all (const string sc, const int ens, const int bin_size, const int time,
                                const string ich_name[3][3], cmatrix Pmat, const string och_name[3]);
   
   void in_out_NBSdata_2x2_spin(const string sc, const int ens, const int bin_size, const int time,
                                const string ich_name[3][3], cmatrix Pmat, const string och_name[3]);
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main() {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   int    bin_size[4] = {0, 57, 40, 45};
   string had_name[4] = {"proton_CG05_CG05", "Lambda_CG05_CG05", "Sigma_CG05_CG05", "Xi_CG05_CG05"};
   string sc[2] = {"s", "c"};
   
   string ich_name_str[2][3][3] =
   {
      {
         {"BBwave.dir.S1.00", "BBwave.dir.S1.01", "BBwave.dir.S1.02"},
         {"BBwave.dir.S1.03", "BBwave.dir.S1.04", "BBwave.dir.S1.05"},
         {"BBwave.dir.S1.06", "BBwave.dir.S1.07", "BBwave.dir.S1.08"}
      },{
         {"BBwave.dir.S2.00", "BBwave.dir.S2.01", "BBwave.dir.S2.02"},
         {"BBwave.dir.S2.06", "BBwave.dir.S2.07", "BBwave.dir.S2.08"},
         {"BBwave.dir.S2.12", "BBwave.dir.S2.13", "BBwave.dir.S2.14"}
      }
   };
   
   string och_name_str[2][2][3] =
   {
      {
         {"LsN12", "SsN12", "SsN32"},
         {"LcN12", "ScN12", "ScN32"}
      },{
         {"LsLs0", "XssN0", "XssN1"},
         {"LcLc0", "XccN0", "XccN1"}
      }
   };
   
   for (int iens=1; iens<=3; iens++) for (int isc=0; isc<2; isc++) {
      for (int ihad=0; ihad<4; ihad++) in_out_CORRdata(sc[isc], iens, bin_size[iens], had_name[ihad]);
      
      for (int itime=max_time; itime>=min_time; itime--) {
         in_out_NBSdata_2x2_all (sc[isc], iens, bin_size[iens], itime, ich_name_str[0], Pmat_LxN , och_name_str[0][isc]);
         in_out_NBSdata_2x2_all (sc[isc], iens, bin_size[iens], itime, ich_name_str[1], Pmat_LxLx, och_name_str[1][isc]);
         in_out_NBSdata_2x2_spin(sc[isc], iens, bin_size[iens], itime, ich_name_str[0], Pmat_LxN , och_name_str[0][isc]);
         in_out_NBSdata_2x2_spin(sc[isc], iens, bin_size[iens], itime, ich_name_str[1], Pmat_LxLx, och_name_str[1][isc]);
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args() {
      // Isospin projection matrix for S=-1 channel
      Pmat_LxN (0,0) = 1.0; Pmat_LxN (0,1) =           0.0 ; Pmat_LxN (0,2) =           0.0 ;
      Pmat_LxN (1,0) = 0.0; Pmat_LxN (1,1) = -sqrt(1.0/3.0); Pmat_LxN (1,2) = +sqrt(2.0/3.0);
      Pmat_LxN (2,0) = 0.0; Pmat_LxN (2,1) = +sqrt(2.0/3.0); Pmat_LxN (2,2) = +sqrt(1.0/3.0);
      
      // Isospin projection matrix for S=-2 channel
      Pmat_LxLx(0,0) = 1.0; Pmat_LxLx(0,1) =           0.0 ; Pmat_LxLx(0,2) =           0.0 ;
      Pmat_LxLx(1,0) = 0.0; Pmat_LxLx(1,1) = -sqrt(1.0/2.0); Pmat_LxLx(1,2) = +sqrt(1.0/2.0);
      Pmat_LxLx(2,0) = 0.0; Pmat_LxLx(2,1) = +sqrt(1.0/2.0); Pmat_LxLx(2,2) = +sqrt(1.0/2.0);
      
      return true;
   }
   ////// ------------------------------------------------------ //////
   void in_out_CORRdata(const string sc, const int ens, const int bin_size, const string hadron) {
      anaHAL::NameList conf(path_conf_list(sc, ens, bin_size));
      STATISTICS<ComplexField_T> corr(conf.Nlist());
      
      printf("Reading Correlator: %s quark, ens=%d [%s]...\n", sc.c_str(), ens, hadron.c_str());
      for (int i=0; i<conf.Nlist(); i++) {
         corr(i).mem_alloc(Tsize);
         corr(i).input_data_corr(path_corr_res(sc, ens, bin_size, conf(i), hadron), true);
      }
      printf("Reading Correlator: %s quark, ens=%d [%s]... end\n", sc.c_str(), ens, hadron.c_str());
      /*
       printf("Make Jack-knife samples... \n");
       corr.make_JK_sample(bin_size);
       printf("Make Jack-knife samples... end\n");
       */
      printf("Writing Correlator: %s quark, ens=%d [%s]...\n", sc.c_str(), ens, hadron.c_str());
      corr.output_data_bin(path_corr_bin(sc, ens, bin_size, hadron));
      printf("Writing Correlator: %s quark, ens=%d [%s]... end\n", sc.c_str(), ens, hadron.c_str());
   }
   ////// ------------------------------------------------------ //////
   void in_out_NBSdata_2x2_all(const string sc, const int ens, const int bin_size, const int time,
                               const string ich_name[3][3], cmatrix Pmat, const string och_name[3]) {
      anaHAL::NameList conf(path_conf_list(sc, ens, bin_size));
      STATISTICS<ComplexField_AXYZB> wave_2x2[2][2];
      STATISTICS<ComplexField_AXYZB> wave_1x1(conf.Nlist());
      for (int i=0; i<2; i++) for (int j=0; j<2; j++) wave_2x2[i][j].mem_alloc(conf.Nlist());
      
      ComplexField_AXYZB wave_in(2*2, Lsize, 2*2);
      ComplexFieldMatrix<ComplexField_AXYZB> wave_mat(3);
      
      printf("Reading NBSwave: %s quark, ens=%d [it=%d]...\n", sc.c_str(), ens, time);
      for (int i=0; i<conf.Nlist(); i++) {
         for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
            wave_in.input_data_comp(path_wave_res(sc, ens, bin_size, ich_name[ich][jch], conf(i), time));
            wave_mat(ich, jch) = wave_in;
         }
         wave_mat = Pmat * wave_mat * Pmat.T();
         
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) wave_2x2[ich][jch](i) = wave_mat(ich,jch);
         wave_1x1(i) = wave_mat(2,2);
         
         printf("Reading NBSwave: %s quark, ens=%d [it=%d]... end: conf=%d\n", sc.c_str(), ens, time, i);
      }
      /*
       printf("Make Jack-knife samples... \n");
       for (int i=0; i<2; i++) for (int j=0; j<2; j++) wave_2x2[i][j].make_JK_sample(bin_size);
       wave_1x1.make_JK_sample(bin_size);
       printf("Make Jack-knife samples... end\n");
       */
      printf("Writing NBSwave: %s quark, ens=%d [it=%d]...\n", sc.c_str(), ens, time);
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         wave_2x2[ich][jch].output_data_bin(path_wave_bin(sc, ens, bin_size, och_name[ich], och_name[jch],
                                                          "allspin", time));
      }
      wave_1x1.output_data_bin(path_wave_bin(sc, ens, bin_size, och_name[2], och_name[2], "allspin", time));
      printf("Writing NBSwave: %s quark, ens=%d [it=%d]... end\n", sc.c_str(), ens, time);
   }
   ////// ------------------------------------------------------ //////
   void in_out_NBSdata_2x2_spin(const string sc, const int ens, const int bin_size, const int time,
                                const string ich_name[3][3], cmatrix Pmat, const string och_name[3]) {
      anaHAL::NameList conf(path_conf_list(sc, ens, bin_size));
      STATISTICS<ComplexField_XYZ> wave_2x2[2][2][2];  // [2]=spin, 0 and 1
      STATISTICS<ComplexField_XYZ> wave_1x1[2];
      for (int ispin=0; ispin<2; ispin++) {
         for (int i=0; i<2; i++) for (int j=0; j<2; j++) wave_2x2[ispin][i][j].mem_alloc(conf.Nlist());
         wave_1x1[ispin].mem_alloc(conf.Nlist());
      }
      ComplexField_AXYZB wave_in(2*2, Lsize, 2*2);
      ComplexFieldMatrix<ComplexField_XYZ> wave_mat[2]; // [2]=spin, 0 and 1
      for (int ispin=0; ispin<2; ispin++)  wave_mat[ispin].init(3);
      
      printf("Reading NBSwave: %s quark, ens=%d [it=%d]...\n", sc.c_str(), ens, time);
      for (int i=0; i<conf.Nlist(); i++) {
         for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
            wave_in.input_data_comp(path_wave_res(sc, ens, bin_size, ich_name[ich][jch], conf(i), time));
            wave_mat[0](ich,jch)  = wave_in.spin_proj(HH_OctOct, 0, 0, HH_OctOct, 0, 0);
            wave_mat[1](ich,jch)  = wave_in.spin_proj(HH_OctOct, 1,+1, HH_OctOct, 1,+1);
            wave_mat[1](ich,jch) += wave_in.spin_proj(HH_OctOct, 1, 0, HH_OctOct, 1, 0);
            wave_mat[1](ich,jch) += wave_in.spin_proj(HH_OctOct, 1,-1, HH_OctOct, 1,-1);
            wave_mat[1](ich,jch) /= 3.0;
         }
         for (int ispin=0; ispin<2; ispin++) {
            wave_mat[ispin] = Pmat * wave_mat[ispin] * Pmat.T();
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               wave_2x2[ispin][ich][jch](i) = wave_mat[ispin](ich,jch);
            }
            wave_1x1[ispin](i) = wave_mat[ispin](2,2);
         }
         printf("Reading NBSwave: %s quark, ens=%d [it=%d]... end: conf=%d\n", sc.c_str(), ens, time, i);
      }
      /*
       printf("Make Jack-knife samples... \n");
       for (int ispin=0; ispin<2; ispin++) {
       for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
       wave_2x2[ispin][ich][jch].make_JK_sample(bin_size);
       }
       wave_1x1[ispin].make_JK_sample(bin_size);
       }
       printf("Make Jack-knife samples... end\n");
       */
      printf("Writing NBSwave: %s quark, ens=%d [it=%d]...\n", sc.c_str(), ens, time);
      for (int ispin=0; ispin<2; ispin++) {
         char spin_c[16]; snprintf(spin_c, sizeof(spin_c), "spin_%d_", ispin); string spin_str(spin_c);
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            wave_2x2[ispin][ich][jch].output_data_bin(path_wave_bin(sc, ens, bin_size, och_name[ich],
                                                                    och_name[jch], spin_str, time));
         }
         wave_1x1[ispin].output_data_bin(path_wave_bin(sc, ens, bin_size, och_name[2], och_name[2], spin_str, time));
      }
      printf("Writing NBSwave: %s quark, ens=%d [it=%d]... end\n", sc.c_str(), ens, time);
   }
}
