#include <ComplexField_Sub.h>
#include <Potential.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   double lattice_spacing = 0.0907;
   
   int min_time = 7;
   int max_time = 15;
   
   string fbase = "/home/miiya369/data/pacs-cs/analysis_20180129";
   
   string path_corr_bin(string sc, int ens, int bin_size, string hadron) {
      char tmp_c[256];
      snprintf(tmp_c, sizeof(tmp_c), "%s_quark/ens%d/bin/correlator.%s.bin_size%02d",
               sc.c_str(), ens, hadron.c_str(), bin_size);
      return (fbase+"/"+tmp_c);
   }
   string path_wave_bin(string sc, int ens, int bin_size, string ch_snk, string ch_src, int spin, int time) {
      char tmp_c[256];
      snprintf(tmp_c, sizeof(tmp_c), "%s_quark/ens%d/bin/NBSwave.%s_%s.spin_%d_.t%02d.bin_size%2d",
               sc.c_str(), ens, ch_snk.c_str(), ch_src.c_str(), spin, time, bin_size);
      return (fbase+"/"+tmp_c);
   }
   string path_data_txt(string dhead, string sc, int ens, int bin_size, string ch_snk, string ch_src,
                        int spin, int time) {
      char tmp_c[256];
      snprintf(tmp_c, sizeof(tmp_c), "%s_quark/ens%d/%s.%s_%s.spin_%d_.t%02d.bin_size%2d.gnu",
               sc.c_str(), ens, dhead.c_str(), ch_snk.c_str(), ch_src.c_str(), spin, time, bin_size);
      return (fbase+"/"+tmp_c);
   }
   string path_data_fit(string dhead, string sc, int ens, int bin_size, string ch_snk, string ch_src,
                        int spin, int time) {
      char tmp_c[256];
      snprintf(tmp_c, sizeof(tmp_c), "%s_quark/ens%d/%s.%s_%s.spin_%d_.t%02d.bin_size%2d.fit",
               sc.c_str(), ens, dhead.c_str(), ch_snk.c_str(), ch_src.c_str(), spin, time, bin_size);
      return (fbase+"/"+tmp_c);
   }
   
   double mass_nucl[3] = {0.7268258335022685, 0.6428135620727219, 0.5586370323829655};
   double mass_Lams[3] = {0.7544946295997861, 0.6861597470193608, 0.6169201873251138};
   double mass_Sigs[3] = {0.7616474016001490, 0.6993680845652021, 0.6413022436222968};
   double mass_Xiss[3] = {0.7854166365174760, 0.7351517471492632, 0.6885019618317600};
   double mass_Lamc[3] = {1.2339722254604193, 1.1743709112922682, 1.1187757095394690};
   double mass_Sigc[3] = {1.2779744472084695, 1.2291281536138188, 1.1834924665483211};
   double mass_Xicc[3] = {1.7475740895660838, 1.7132121754793537, 1.6888517284982165};
   
   double sqZP_nucl[3] = {0.0638380021363467, 0.0479060440234551, 0.0374541258735672};
   double sqZP_Lams[3] = {0.0680915259596637, 0.0537393111098323, 0.0423508918810769};
   double sqZP_Sigs[3] = {0.0672010952825638, 0.0534853284888701, 0.0425273646365216};
   double sqZP_Xiss[3] = {0.0720090097129906, 0.0595294964305716, 0.0495200977265162};
   double sqZP_Lamc[3] = {0.1176169411797812, 0.0932302955882113, 0.0778753048612198};
   double sqZP_Sigc[3] = {0.1025768292035812, 0.0851941294080675, 0.0685015596460554};
   double sqZP_Xicc[3] = {0.1942732211319675, 0.1728605502306786, 0.1629757104586813};
   
   STATISTICS<ComplexField_T  > corr[2][2];
   STATISTICS<ComplexField_XYZ> wave[3][2][2];
   STATISTICS<ComplexField_XYZ>  pot[2][2];
   
   bool set_args();
   
   void input_corrdata_from_bin(const string sc, const int ens, const int bin_size,
                                const string had_name[2][2]);
   
   void input_NBSdata_from_bin (const string sc, const int ens, const int bin_size, const int time,
                                const int spin, const string ich_name[2]);
   
   void potential_calculation  (const int time, const double mass[2][2], const double sqZP[2][2]);
   
   void output_XYZdata_to_text (const string sc, const int ens, const int bin_size, const int time,
                                const int spin, const string ich_name[2], const string dir_name);
   void output_XYZdata_to_bin  (const string sc, const int ens, const int bin_size, const int time,
                                const int spin, const string ich_name[2], const string dir_name);
   
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main() {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   int    bin_size[4] = {0, 57, 40, 45};
   string sc[2] = {"s", "c"};
   
   string hadron_names_str[2][2][2] =
   {
      {
         {"proton_CG05_CG05", "Lambda_CG05_CG05"},
         {"proton_CG05_CG05",  "Sigma_CG05_CG05"}
      },{
         {"Lambda_CG05_CG05", "Lambda_CG05_CG05"},
         {"proton_CG05_CG05",     "Xi_CG05_CG05"}
      }
   };
   
   string ich_name_str[2][2][2] =
   {
      {
         {"LsN12", "SsN12"},
         {"LcN12", "ScN12"}
      },{
         {"LsLs0", "XssN0"},
         {"LcLc0", "XccN0"}
      }
   };
   
   string odir_name_str[2] = {"LxN", "LxLx"};
   
   for (int iens=1; iens<=3; iens++) {
      double hadron_masses[2][2][2][2] =
      {
         {
            {
               {mass_nucl[iens-1], mass_Lams[iens-1]},
               {mass_nucl[iens-1], mass_Sigs[iens-1]}
            },{
               {mass_nucl[iens-1], mass_Lamc[iens-1]},
               {mass_nucl[iens-1], mass_Sigc[iens-1]}
            }
         },{
            {
               {mass_Lams[iens-1], mass_Lams[iens-1]},
               {mass_nucl[iens-1], mass_Xiss[iens-1]}
            },{
               {mass_Lamc[iens-1], mass_Lamc[iens-1]},
               {mass_nucl[iens-1], mass_Xicc[iens-1]}
            }
         }
      };
      
      double hadron_sqZP[2][2][2][2] =
      {
         {
            {
               {sqZP_nucl[iens-1], sqZP_Lams[iens-1]},
               {sqZP_nucl[iens-1], sqZP_Sigs[iens-1]}
            },{
               {sqZP_nucl[iens-1], sqZP_Lamc[iens-1]},
               {sqZP_nucl[iens-1], sqZP_Sigc[iens-1]}
            }
         },{
            {
               {sqZP_Lams[iens-1], sqZP_Lams[iens-1]},
               {sqZP_nucl[iens-1], sqZP_Xiss[iens-1]}
            },{
               {sqZP_Lamc[iens-1], sqZP_Lamc[iens-1]},
               {sqZP_nucl[iens-1], sqZP_Xicc[iens-1]}
            }
         }
      };
      for (int isc=0; isc<2; isc++) for (int ich=0; ich<2; ich++) {
         input_corrdata_from_bin(sc[isc], iens, bin_size[iens], hadron_names_str[ich]);
         
         for (int ispin=0; ispin<=1; ispin++) for (int itime=max_time; itime>=min_time; itime--) {
            input_NBSdata_from_bin(sc[isc], iens, bin_size[iens], itime, ispin, ich_name_str[ich][isc]);
            potential_calculation (itime, hadron_masses[ich][isc], hadron_sqZP[ich][isc]);
            
            output_XYZdata_to_text(sc[isc], iens, bin_size[iens], itime, ispin,
                                   ich_name_str[ich][isc], odir_name_str[ich]);
            output_XYZdata_to_bin (sc[isc], iens, bin_size[iens], itime, ispin,
                                   ich_name_str[ich][isc], odir_name_str[ich]);
         }
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args() {
      
      return true;
   }
   ////// ------------------------------------------------------ //////
   void input_corrdata_from_bin(const string sc, const int ens, const int bin_size,
                                const string had_name[2][2]) {
      for (int ich=0; ich<2; ich++)
         for (int jch=0; jch<2; jch++)
            corr[ich][jch].input_data_bin(path_corr_bin(sc, ens, bin_size, had_name[ich][jch]));
   }
   ////// ------------------------------------------------------ //////
   void input_NBSdata_from_bin (const string sc, const int ens, const int bin_size, const int time,
                                const int spin, const string ich_name[2]) {
      for (int dt=-1; dt<=1; dt++)
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++)
            wave[dt+1][ich][jch].input_data_bin(path_wave_bin(sc, ens, bin_size, ich_name[ich],
                                                              ich_name[jch], spin, time+dt));
   }
   ////// ------------------------------------------------------ //////
   void potential_calculation(const int time, const double mass[2][2], const double sqZP[2][2]) {
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++)
         pot[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
      
#pragma omp parallel for
      for (int i=0; i<wave[0][0][0].Ndata(); i++) {
         ComplexField_XYZ tmp_pot[2*2], tmp_Rcor[3][2*2];
         cdouble          tmp_fac[2*2];
         double           mass1[2], mass2[2];
         
         for (int ich=0; ich<2; ich++) {
            mass1[ich] = mass[ich][0];
            mass2[ich] = mass[ich][1];
            
            for (int jch=0; jch<2; jch++) {
               int ijch = jch + 2 * ich;
               tmp_fac[ijch] =
                  ((sqZP[ich][0] * sqZP[ich][1] * corr[jch][0](i)(time) * corr[jch][1](i)(time)) /
                   (sqZP[jch][0] * sqZP[jch][1] * corr[ich][0](i)(time) * corr[ich][1](i)(time)));
               
               for (int dt=-1; dt<=+1; dt++)
                  tmp_Rcor[dt+1][ijch] = (wave[dt+1][ich][jch](i).rot_proj(ROT_REP_A1) /
                                          (corr[ich][0](i)(time+dt) * corr[ich][1](i)(time+dt)));
            }
         }
         Potential::calc_CCpotential_T2(tmp_pot, tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                        mass1, mass2, tmp_fac, 2);
         
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            int ijch = jch + 2 * ich;
            pot[ich][jch](i) = tmp_pot[ijch];
         }
      }
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_text(const string sc, const int ens, const int bin_size, const int time,
                               const int spin, const string ich_name[2], const string dir_name) {
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         wave[1][ich][jch].output_data_err(path_data_txt(dir_name+"/wave/NBS", sc, ens, bin_size,
                                                         ich_name[ich], ich_name[jch],
                                                         spin, time), lattice_spacing, true);
         pot [ich][jch].output_data_err(path_data_txt(dir_name+"/pot/Pot"  , sc, ens, bin_size,
                                                      ich_name[ich], ich_name[jch],
                                                      spin, time), lattice_spacing, true);
      }
   }
   ////// ------------------------------------------------------ //////
   void output_XYZdata_to_bin(const string sc, const int ens, const int bin_size, const int time,
                              const int spin, const string ich_name[2], const string dir_name) {
      for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
         pot[ich][jch].output_data_bin_reduce(path_data_fit(dir_name+"/pot/Pot", sc, ens, bin_size,
                                                            ich_name[ich], ich_name[jch],
                                                            spin, time), lattice_spacing, true);
      }
   }
}
