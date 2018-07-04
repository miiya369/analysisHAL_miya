#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int    Tsize = 48;
   int    Lsize = 48;
   
   int min_time = 5;
   int max_time = 18;
   
   int bin_size = 40;
   int Nconf    = 800      / 4; // 4 = x,y,z,t-rotation
   int tmp_div  = bin_size / 4; // these shold be dividable
   
   string fbase = "/home/miiya369/data/yku/48x48/analysis_20180521/LN";
   
   string conf_str(int iconf, int ixyzt) {
      string cbase   = "RC48x48_B1900Kud01373316Ks01367526C1715-gM-GFXD_PREC1.e-14";
      string xyzt[4] = {"tyzx", "xtzy", "xytz", "xyzt"};
      int c_start    =  41;
      
      char conf_c[256];
      snprintf(conf_c, sizeof(conf_c), "%s_%05d0.rot_%s", cbase.c_str(), iconf+c_start, xyzt[ixyzt].c_str());
      return string(conf_c);
   }
   string path_corr_res(string ibase, string conf, string hadron, string src_shift_str) {
      return (ibase+"/correlator.PS.dir/"+conf+"/"+hadron+"_correlator.+"+src_shift_str+"."+conf);
   }
   string path_wave_res(string ibase, string dir_channel, string conf, int time, string src_shift_str) {
      char time_c[8]; snprintf(time_c, sizeof(time_c), "+%03d", time); string time_str(time_c);
      return (ibase+"/"+dir_channel+"/"+conf+"/NBSwave."+time_str+"+"+src_shift_str+"."+conf+".NUC_CG05.NUC_CG05");
   }
   string path_corr_bin(string obase, string we, int bsize, string hadron) {
      char bsize_c[8]; snprintf(bsize_c, sizeof(bsize_c), "%02d", bsize); string bsize_str(bsize_c);
      return (obase+"/"+we+"/bin/correlator."+we+"."+hadron+".bin_size"+bsize_str);
   }
   string path_wave_bin(string obase, string we, int bsize, string ch_snk, string ch_src, string spin, int time) {
      char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",  time); string  time_str( time_c);
      char bsize_c[8]; snprintf(bsize_c, sizeof(bsize_c), "%02d", bsize); string bsize_str(bsize_c);
      return (obase+"/"+we+"/bin/NBSwave."+we+"."+ch_snk+"_"+ch_src+"."+spin+".t"+time_str+".bin_size"+bsize_str);
   }
   
   string wall_expo_str[2] = {"wall", "expo"};
   string src_shift_str[2] = {"B00.000.000.000", "A00.006_018_A00.A00.A00"};
   
   string had_name_str[3] = {"proton_CG05_CG05", "Lambda_CG05_CG05", "Sigma_CG05_CG05"};
   string ich_name_str[3][3] =
   {
      {"BBwave.dir.S1.00", "BBwave.dir.S1.01", "BBwave.dir.S1.02"},
      {"BBwave.dir.S1.03", "BBwave.dir.S1.04", "BBwave.dir.S1.05"},
      {"BBwave.dir.S1.06", "BBwave.dir.S1.07", "BBwave.dir.S1.08"}
   };
   string och_name_str[3] = {"LamN12", "SigN12", "SigN32"};
   
   cmatrix Pmat_LN_SN(3);
   
   bool set_args();
   void io_corr       (const string we, const string src_shift, const string hadron);
   void io_NBS_2x2_all(const string we, const string src_shift, const int time,
                       const string ich_name[3][3], cmatrix Pmat, const string och_name[3]);
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main() {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   for (int iwe=0; iwe<2; iwe++) {
      for (int ihad=0; ihad<3; ihad++) io_corr(wall_expo_str[iwe], src_shift_str[iwe], had_name_str[ihad]);
      
      for (int itime=max_time; itime>=min_time; itime--)
         io_NBS_2x2_all(wall_expo_str[iwe], src_shift_str[iwe], itime, ich_name_str, Pmat_LN_SN, och_name_str);
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args() {
      // Isospin projection matrix for S=-1 channel
      Pmat_LN_SN(0,0) = 1.0; Pmat_LN_SN(0,1) =           0.0 ; Pmat_LN_SN(0,2) =           0.0 ;
      Pmat_LN_SN(1,0) = 0.0; Pmat_LN_SN(1,1) = -sqrt(1.0/3.0); Pmat_LN_SN(1,2) = +sqrt(2.0/3.0);
      Pmat_LN_SN(2,0) = 0.0; Pmat_LN_SN(2,1) = +sqrt(2.0/3.0); Pmat_LN_SN(2,2) = +sqrt(1.0/3.0);
      
      return true;
   }
   ////// ------------------------------------------------------ //////
   void io_corr(const string we, const string src_shift, const string hadron) {
      STATISTICS<ComplexField_T> corr(4*Nconf); // 4 = x,y,z,t-rotation
      
      printf("Reading Correlator: src=%s [%s]... ", we.c_str(), hadron.c_str()); fflush(stdout);
      for (int iconf=0; iconf<Nconf; iconf++) for (int ixyzt=0; ixyzt<4; ixyzt++) {
         int i = ixyzt + 4 * iconf;
         corr(i).mem_alloc(Tsize);
         corr(i).input_data_corr(path_corr_res(fbase+"/"+we+"/data", conf_str(iconf,ixyzt), hadron, src_shift), true);
      }
      printf("end\n");
      
      printf("Making  JK-samples: src=%s [%s]... ", we.c_str(), hadron.c_str()); fflush(stdout);
      corr.make_JK_sample(bin_size);
      printf("end\n");
      
      printf("Writing Correlator: src=%s [%s]... ", we.c_str(), hadron.c_str()); fflush(stdout);
      corr.output_data_bin(path_corr_bin(fbase, we, bin_size, hadron));
      printf("end\n");
   }
   ////// ------------------------------------------------------ //////
   void io_NBS_2x2_all(const string we, const string src_shift, const int time,
                       const string ich_name[3][3], cmatrix Pmat, const string och_name[3]) {
      STATISTICS<ComplexField_AXYZB> wave_2x2[2][2];
      STATISTICS<ComplexField_AXYZB> wave_1x1(Nconf/tmp_div);
      for (int i=0; i<2; i++) for (int j=0; j<2; j++) wave_2x2[i][j].mem_alloc(Nconf/tmp_div);
      
      ComplexField_AXYZB wave_in(2*2, Lsize, 2*2);
      ComplexFieldMatrix<ComplexField_AXYZB> wave_mat(3);
      for (int ijch=0; ijch<3*3; ijch++) wave_mat(ijch).mem_alloc(2*2, Lsize, 2*2);
      
      for (int iconf=0; iconf<Nconf/tmp_div; iconf++) {
         for (int ich=0; ich<3; ich++) for (int jch=0; jch<3; jch++) {
            wave_mat(ich, jch) = 0.0;
            for (int jconf=0; jconf<tmp_div; jconf++) for (int ixyzt=0; ixyzt<4; ixyzt++) {
               int i = jconf + tmp_div * iconf;
               wave_in.input_data_comp(path_wave_res(fbase+"/"+we+"/data", ich_name[ich][jch], conf_str(i,ixyzt), time, src_shift));
               wave_mat(ich, jch) += wave_in;
            }
            wave_mat(ich, jch) /= (double)bin_size;
         }
         wave_mat = Pmat * wave_mat * Pmat.T();
         
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) wave_2x2[ich][jch](iconf) = wave_mat(ich,jch);
         wave_1x1(iconf) = wave_mat(2,2);
         printf("Reading NBSwave: src=%s [it=%d, conf=%d]... end\n", we.c_str(), time, iconf);
      }
      
      printf("Making JKsample: src=%s [it=%d]... ", we.c_str(), time); fflush(stdout);
      for (int i=0; i<2; i++) for (int j=0; j<2; j++) wave_2x2[i][j].make_JK_sample(1);
      wave_1x1.make_JK_sample(1);
      printf("end\n");
      
      printf("Writing NBSwave: src=%s [it=%d]... ", we.c_str(), time); fflush(stdout);
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         wave_2x2[ich][jch].output_data_bin(path_wave_bin(fbase, we, bin_size, och_name[ich], och_name[jch], "allspin", time));
      }
      wave_1x1.output_data_bin(path_wave_bin(fbase, we, bin_size, och_name[2], och_name[2], "allspin", time));
      printf("end\n");
   }
}
