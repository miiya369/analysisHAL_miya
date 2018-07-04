#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>
#include <MatrixFunc.h>
#include <Potential.h>

///////// =========================== MAIN PART =========================== /////////
int main() {
   time_t stime, etime; time(&stime);
   
   /// Set parameterschar ///
   double lattice_spacing = 0.0907;
   int    bin_size        = 57;
   int    min_time        = 7;
   int    max_time        = 15;
   
   string quark = "s";
   
   string ibase = "/home/miiya369/data/pacs-cs/analysis_20180515/"+quark+"_quark/ens1/bin";
   string obase = "/home/miiya369/data/pacs-cs/analysis_20180515/"+quark+"_quark/ens1/LxN/pot";
   
   double mass_nucl  = 0.726840;
   
   double mass_Lam   = 0.754507;
   double mass_Sig   = 0.761660;
   /*
   double mass_Lam   = 1.233977;
   double mass_Sig   = 1.277982;
   */
   double sq_ZP_nucl  = 0.0638732520686627;
   
   double sq_ZP_Lam   = 0.0681264891598947;
   double sq_ZP_Sig   = 0.0672251624030046;
   /*
   double sq_ZP_Lam   = 0.1176390522262307;
   double sq_ZP_Sig   = 0.1025945201743403;
   */
   
   char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
   
   /////////////////////////////////////////////////////////////////
   /// For I=1/2 , J=0 coupled channel (e.g. Lambda N - Sigma N) ///
   {
      printf("@@@ I=1/2 , J=0 coupled channel calculating... "); fflush(stdout);
      string ihad_str[2][2] = {
         {"proton_CG05_CG05", "proton_CG05_CG05"},
         {"Lambda_CG05_CG05", "Sigma_CG05_CG05"},
      };
      double mass[2][2] = {
         {mass_nucl, mass_nucl},
         {mass_Lam,  mass_Sig}
      };
      cdouble ZP[2] = {sq_ZP_nucl*sq_ZP_Lam, sq_ZP_nucl*sq_ZP_Sig};
      
      string ich_str[2] = {"L"+quark+"N12_", "S"+quark+"N12_"};
      
      STATISTICS<ComplexField_T  > corr[2][2];
      STATISTICS<ComplexField_XYZ> wave[3][2][2];
      STATISTICS<ComplexField_XYZ>  pot[2][2];
      STATISTICS<ComplexField_XYZ>  pot_odiag;
      
      for (   int ihad=0; ihad<2; ihad++)
         for (int  ich=0;  ich<2; ich++ )
            corr[ihad][ich].input_data_bin(ibase+"/correlator."+ihad_str[ihad][ich]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            for (   int ich=0; ich<2; ich++)
               for (int jch=0; jch<2; jch++)
                  wave[dt+1][ich][jch].input_data_bin(ibase+"/NBSwave."+ich_str[ich]+"_"+ich_str[jch]+
                                                      ".spin_0_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         for (   int ich=0; ich<2; ich++)
            for (int jch=0; jch<2; jch++)
               pot[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
         
         pot_odiag.mem_alloc(pot[0][0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0][0][0].Ndata(); i++) {
            ComplexField_XYZ tmp_pot[2*2], tmp_Rcor[3][2*2];
            cdouble          tmp_fac[2*2];
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int  ijch = jch + 2*ich;
               tmp_fac[ijch] = ((ZP[ich] * corr[0][jch](i)(itime) * corr[1][jch](i)(itime)) /
                                (ZP[jch] * corr[0][ich](i)(itime) * corr[1][ich](i)(itime)));
               
               for (int dt=-1; dt<=+1; dt++)
                  tmp_Rcor[dt+1][ijch] = (wave[dt+1][ich][jch](i).rot_proj(ROT_REP_A1) /
                                          (corr[0][ich](i)(itime+dt) * corr[1][ich](i)(itime+dt)));
            }
            Potential::calc_CCpotential_T2(tmp_pot, tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                           mass[0], mass[1], tmp_fac, 2);
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int  ijch = jch + 2*ich;
               pot[ich][jch](i) = tmp_pot[ijch];
            }
            pot_odiag(i) = (pot[0][1](i) + pot[1][0](i)) / 2.0;
         }
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            string ofname_base = (obase+"/Pot."+ich_str[ich]+"_"+ich_str[jch]+
                                  ".V_C.J_0.t"+time_str+".bin_size"+Bsize_str);
            
            pot[ich][jch].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
            pot[ich][jch].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
            pot[ich][jch].output_data_bin       (ofname_base+".bin");
         }
         string ofname_base = obase+"/Pot."+quark+"_offdiag_ave.V_C.J_0.t"+time_str+".bin_size"+Bsize_str;
         
         pot_odiag.output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
         pot_odiag.output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
      }
      printf("END\n");
   }
   
   //////////////////////////////////////////////////////
   /// For I=1/2 , J=0 single channel (e.g. Lambda N) ///
   {
      printf("@@@ I=1/2 , J=0 single  channel calculating... "); fflush(stdout);
      
      string ihad_str[2] = {"proton_CG05_CG05", "Lambda_CG05_CG05"};
      double     mass[2] = {         mass_nucl,           mass_Lam};
      
      string ich_str = "L"+quark+"N12_";
      
      STATISTICS<ComplexField_T  > corr[2];
      STATISTICS<ComplexField_XYZ> wave[3];
      STATISTICS<ComplexField_XYZ>  pot;
      
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(ibase+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            wave[dt+1].input_data_bin(ibase+"/NBSwave."+ich_str+"_"+ich_str+
                                      ".spin_0_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         pot.mem_alloc(wave[0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0].Ndata(); i++) {
            ComplexField_XYZ tmp_Rcor[3];
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1] = (wave[dt+1](i).rot_proj(ROT_REP_A1) /
                                 (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
            
            pot(i) = Potential::get_potential_T2(tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1]);
         }
         string ofname_base = obase+"/Pot."+ich_str+"_single.V_C.J_0.t"+time_str+".bin_size"+Bsize_str;
         
         pot.output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
         pot.output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
         pot.output_data_bin       (ofname_base+".bin");
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////
   /// For I=3/2 , J=0 single channel (e.g. Sigma N) ///
   {
      printf("@@@ I=3/2 , J=0 single  channel calculating... "); fflush(stdout);
      
      string ihad_str[2] = {"proton_CG05_CG05", "Sigma_CG05_CG05"};
      double     mass[2] = {         mass_nucl,          mass_Sig};
      
      string ich_str = "S"+quark+"N32_";
      
      STATISTICS<ComplexField_T  > corr[2];
      STATISTICS<ComplexField_XYZ> wave[3];
      STATISTICS<ComplexField_XYZ>  pot;
      
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(ibase+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            wave[dt+1].input_data_bin(ibase+"/NBSwave."+ich_str+"_"+ich_str+
                                      ".spin_0_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         pot.mem_alloc(wave[0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0].Ndata(); i++) {
            ComplexField_XYZ tmp_Rcor[3];
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1] = (wave[dt+1](i).rot_proj(ROT_REP_A1) /
                                 (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
            
            pot(i) = Potential::get_potential_T2(tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1]);
         }
         string ofname_base = obase+"/Pot."+ich_str+"_single.V_C.J_0.t"+time_str+".bin_size"+Bsize_str;
         
         pot.output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
         pot.output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
         pot.output_data_bin       (ofname_base+".bin");
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////////////////
   /// For I=1/2 , J=1 coupled channel (e.g. Lambda N - Sigma N) ///
   {
      printf("@@@ I=1/2 , J=1 coupled channel calculating... "); fflush(stdout);
      
      string ihad_str[2][2] = {
         {"proton_CG05_CG05", "proton_CG05_CG05"},
         {"Lambda_CG05_CG05", "Sigma_CG05_CG05"},
      };
      double mass[2][2] = {
         {mass_nucl, mass_nucl},
         {mass_Lam,  mass_Sig}
      };
      cdouble ZP[2] = {sq_ZP_nucl*sq_ZP_Lam, sq_ZP_nucl*sq_ZP_Sig};
      
      string ich_str[2] = {"L"+quark+"N12_", "S"+quark+"N12_"};
      
      STATISTICS<ComplexField_T    > corr[2][2];
      STATISTICS<ComplexField_AXYZB> wave[3][2][2];
      STATISTICS<ComplexField_XYZ  > potC[2][2], potT[2][2];
      STATISTICS<ComplexField_XYZ>   potC_odiag, potT_odiag;
      
      for (   int ihad=0; ihad<2; ihad++)
         for (int  ich=0;  ich<2; ich++ )
            corr[ihad][ich].input_data_bin(ibase+"/correlator."+ihad_str[ihad][ich]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            for (   int ich=0; ich<2; ich++)
               for (int jch=0; jch<2; jch++)
                  wave[dt+1][ich][jch].input_data_bin(ibase+"/NBSwave."+ich_str[ich]+"_"+ich_str[jch]+
                                                      ".allspin.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         for (   int ich=0; ich<2; ich++)
            for (int jch=0; jch<2; jch++) {
               potC[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
               potT[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
            }
         potC_odiag.mem_alloc(potC[0][0].Ndata());
         potT_odiag.mem_alloc(potT[0][0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0][0][0].Ndata(); i++) {
            ComplexField_AXYZB tmp_Rcor[3][2*2];
            ComplexField_XYZ   tmp_potC[2*2], tmp_potT[2*2];
            cdouble            tmp_fac[2*2];
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int  ijch = jch + 2*ich;
               tmp_fac[ijch] = ((ZP[ich] * corr[0][jch](i)(itime) * corr[1][jch](i)(itime)) /
                                (ZP[jch] * corr[0][ich](i)(itime) * corr[1][ich](i)(itime)));
               
               for (int dt=-1; dt<=+1; dt++)
                  tmp_Rcor[dt+1][ijch] = (wave[dt+1][ich][jch](i) /
                                          (corr[0][ich](i)(itime+dt) * corr[1][ich](i)(itime+dt)));
            }
            Potential::calc_tensor_CCpot_T2(tmp_potC, tmp_potT, tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                            mass[0], mass[1], tmp_fac, 2);
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int  ijch = jch + 2*ich;
               potC[ich][jch](i) = tmp_potC[ijch];
               potT[ich][jch](i) = tmp_potT[ijch];
            }
            potC_odiag(i) = (potC[0][1](i) + potC[1][0](i)) / 2.0;
            potT_odiag(i) = (potT[0][1](i) + potT[1][0](i)) / 2.0;
         }
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            string ofname_base_C = (obase+"/Pot."+ich_str[ich]+"_"+ich_str[jch]+
                                    ".V_C.J_1.t"+time_str+".bin_size"+Bsize_str);
            string ofname_base_T = (obase+"/Pot."+ich_str[ich]+"_"+ich_str[jch]+
                                    ".V_T.J_1.t"+time_str+".bin_size"+Bsize_str);
            
            potC[ich][jch].output_data_err       (ofname_base_C+".gnu", lattice_spacing,  true);
            potC[ich][jch].output_data_bin_reduce(ofname_base_C+".fit", lattice_spacing, false);
            potC[ich][jch].output_data_bin       (ofname_base_C+".bin");
            
            potT[ich][jch].output_data_err       (ofname_base_T+".gnu", lattice_spacing,  true);
            potT[ich][jch].output_data_bin_reduce(ofname_base_T+".fit", lattice_spacing, false);
            potT[ich][jch].output_data_bin       (ofname_base_T+".bin");
         }
         string ofname_base_C = (obase+"/Pot."+quark+"_offdiag_ave.V_C.J_1.t"+
                                 time_str+".bin_size"+Bsize_str);
         string ofname_base_T = (obase+"/Pot."+quark+"_offdiag_ave.V_T.J_1.t"+
                                 time_str+".bin_size"+Bsize_str);
         
         potC_odiag.output_data_err       (ofname_base_C+".gnu", lattice_spacing,  true);
         potC_odiag.output_data_bin_reduce(ofname_base_C+".fit", lattice_spacing, false);
         
         potT_odiag.output_data_err       (ofname_base_T+".gnu", lattice_spacing,  true);
         potT_odiag.output_data_bin_reduce(ofname_base_T+".fit", lattice_spacing, false);
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////
   /// For I=1/2 , J=1 single channel (e.g. Lambda N) ///
   {
      printf("@@@ I=1/2 , J=1 single  channel calculating... "); fflush(stdout);
      
      string ihad_str[2] = {"proton_CG05_CG05", "Lambda_CG05_CG05"};
      double     mass[2] = {         mass_nucl,           mass_Lam};
      
      string ich_str = "L"+quark+"N12_";
      
      STATISTICS<ComplexField_T    > corr[2];
      STATISTICS<ComplexField_AXYZB> wave[3];
      STATISTICS<ComplexField_XYZ  > potC, potT;
      
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(ibase+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            wave[dt+1].input_data_bin(ibase+"/NBSwave."+ich_str+"_"+ich_str+
                                      ".allspin.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         potC.mem_alloc(wave[0].Ndata());
         potT.mem_alloc(wave[0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0].Ndata(); i++) {
            ComplexField_AXYZB tmp_Rcor[3];
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1] = (wave[dt+1](i) / (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
            
            Potential::calc_tensor_pot_T2(potC(i), potT(i), tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1]);
         }
         string ofname_base_C = obase+"/Pot."+ich_str+"_single.V_C.J_1.t"+time_str+".bin_size"+Bsize_str;
         string ofname_base_T = obase+"/Pot."+ich_str+"_single.V_T.J_1.t"+time_str+".bin_size"+Bsize_str;
         
         potC.output_data_err       (ofname_base_C+".gnu", lattice_spacing,  true);
         potC.output_data_bin_reduce(ofname_base_C+".fit", lattice_spacing, false);
         potC.output_data_bin       (ofname_base_C+".bin");
         
         potT.output_data_err       (ofname_base_T+".gnu", lattice_spacing,  true);
         potT.output_data_bin_reduce(ofname_base_T+".fit", lattice_spacing, false);
         potT.output_data_bin       (ofname_base_T+".bin");
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////
   /// For I=3/2 , J=1 single channel (e.g. Sigma N) ///
   {
      printf("@@@ I=3/2 , J=1 single  channel calculating... "); fflush(stdout);
      
      string ihad_str[2] = {"proton_CG05_CG05", "Sigma_CG05_CG05"};
      double     mass[2] = {         mass_nucl,          mass_Sig};
      
      string ich_str = "S"+quark+"N32_";
      
      STATISTICS<ComplexField_T    > corr[2];
      STATISTICS<ComplexField_AXYZB> wave[3];
      STATISTICS<ComplexField_XYZ  > potC, potT;
      
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(ibase+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            wave[dt+1].input_data_bin(ibase+"/NBSwave."+ich_str+"_"+ich_str+
                                      ".allspin.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         potC.mem_alloc(wave[0].Ndata());
         potT.mem_alloc(wave[0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0].Ndata(); i++) {
            ComplexField_AXYZB tmp_Rcor[3];
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1] = (wave[dt+1](i) / (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
            
            Potential::calc_tensor_pot_T2(potC(i), potT(i), tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1]);
         }
         string ofname_base_C = obase+"/Pot."+ich_str+"_single.V_C.J_1.t"+time_str+".bin_size"+Bsize_str;
         string ofname_base_T = obase+"/Pot."+ich_str+"_single.V_T.J_1.t"+time_str+".bin_size"+Bsize_str;
         
         potC.output_data_err       (ofname_base_C+".gnu", lattice_spacing,  true);
         potC.output_data_bin_reduce(ofname_base_C+".fit", lattice_spacing, false);
         potC.output_data_bin       (ofname_base_C+".bin");
         
         potT.output_data_err       (ofname_base_T+".gnu", lattice_spacing,  true);
         potT.output_data_bin_reduce(ofname_base_T+".fit", lattice_spacing, false);
         potT.output_data_bin       (ofname_base_T+".bin");
      }
      printf("END\n");
   }
   
   ////////////////////////////////////////////////////////////////////////////////////
   /// For I=1/2 , J=1 <Effective S-wave> coupled channel (e.g. Lambda N - Sigma N) ///
   {
      printf("@@@ I=1/2 , J=1 <Effective S-wave> coupled channel calculating... "); fflush(stdout);
      string ihad_str[2][2] = {
         {"proton_CG05_CG05", "proton_CG05_CG05"},
         {"Lambda_CG05_CG05", "Sigma_CG05_CG05"},
      };
      double mass[2][2] = {
         {mass_nucl, mass_nucl},
         {mass_Lam,  mass_Sig}
      };
      cdouble ZP[2] = {sq_ZP_nucl*sq_ZP_Lam, sq_ZP_nucl*sq_ZP_Sig};
      
      string ich_str[2] = {"L"+quark+"N12_", "S"+quark+"N12_"};
      
      STATISTICS<ComplexField_T  > corr[2][2];
      STATISTICS<ComplexField_XYZ> wave[3][2][2];
      STATISTICS<ComplexField_XYZ>  pot[2][2];
      STATISTICS<ComplexField_XYZ>  pot_odiag;
      
      for (   int ihad=0; ihad<2; ihad++)
         for (int  ich=0;  ich<2; ich++ )
            corr[ihad][ich].input_data_bin(ibase+"/correlator."+ihad_str[ihad][ich]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            for (   int ich=0; ich<2; ich++)
               for (int jch=0; jch<2; jch++)
                  wave[dt+1][ich][jch].input_data_bin(ibase+"/NBSwave."+ich_str[ich]+"_"+ich_str[jch]+
                                                      ".spin_1_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         for (   int ich=0; ich<2; ich++)
            for (int jch=0; jch<2; jch++)
               pot[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
         
         pot_odiag.mem_alloc(pot[0][0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0][0][0].Ndata(); i++) {
            ComplexField_XYZ tmp_pot[2*2], tmp_Rcor[3][2*2];
            cdouble          tmp_fac[2*2];
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int  ijch = jch + 2*ich;
               tmp_fac[ijch] = ((ZP[ich] * corr[0][jch](i)(itime) * corr[1][jch](i)(itime)) /
                                (ZP[jch] * corr[0][ich](i)(itime) * corr[1][ich](i)(itime)));
               
               for (int dt=-1; dt<=+1; dt++)
                  tmp_Rcor[dt+1][ijch] = (wave[dt+1][ich][jch](i).rot_proj(ROT_REP_A1) /
                                          (corr[0][ich](i)(itime+dt) * corr[1][ich](i)(itime+dt)));
            }
            Potential::calc_CCpotential_T2(tmp_pot, tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                           mass[0], mass[1], tmp_fac, 2);
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               int  ijch = jch + 2*ich;
               pot[ich][jch](i) = tmp_pot[ijch];
            }
            pot_odiag(i) = (pot[0][1](i) + pot[1][0](i)) / 2.0;
         }
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            string ofname_base = (obase+"/Pot."+ich_str[ich]+"_"+ich_str[jch]+
                                  ".eff.3S1.t"+time_str+".bin_size"+Bsize_str);
            
            pot[ich][jch].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
            pot[ich][jch].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
            pot[ich][jch].output_data_bin       (ofname_base+".bin");
         }
         string ofname_base = obase+"/Pot."+quark+"_offdiag_ave.eff.3S1.t"+time_str+".bin_size"+Bsize_str;
         
         pot_odiag.output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
         pot_odiag.output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////////////////////////
   /// For I=1/2 , J=1 <Effective S-wave> single channel (e.g. Lambda N) ///
   {
      printf("@@@ I=1/2 , J=1 <Effective S-wave> single  channel calculating... "); fflush(stdout);
      
      string ihad_str[2] = {"proton_CG05_CG05", "Lambda_CG05_CG05"};
      double     mass[2] = {         mass_nucl,           mass_Lam};
      
      string ich_str = "L"+quark+"N12_";
      
      STATISTICS<ComplexField_T  > corr[2];
      STATISTICS<ComplexField_XYZ> wave[3];
      STATISTICS<ComplexField_XYZ>  pot;
      
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(ibase+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            wave[dt+1].input_data_bin(ibase+"/NBSwave."+ich_str+"_"+ich_str+
                                      ".spin_1_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         pot.mem_alloc(wave[0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0].Ndata(); i++) {
            ComplexField_XYZ tmp_Rcor[3];
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1] = (wave[dt+1](i).rot_proj(ROT_REP_A1) /
                                 (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
            
            pot(i) = Potential::get_potential_T2(tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1]);
         }
         string ofname_base = obase+"/Pot."+ich_str+"_single.eff.3S1.t"+time_str+".bin_size"+Bsize_str;
         
         pot.output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
         pot.output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
         pot.output_data_bin       (ofname_base+".bin");
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////
   /// For I=3/2 , J=0 <Effective S-wave> single channel (e.g. Sigma N) ///
   {
      printf("@@@ I=3/2 , J=0 <Effective S-wave> single  channel calculating... "); fflush(stdout);
      
      string ihad_str[2] = {"proton_CG05_CG05", "Sigma_CG05_CG05"};
      double     mass[2] = {         mass_nucl,          mass_Sig};
      
      string ich_str = "S"+quark+"N32_";
      
      STATISTICS<ComplexField_T  > corr[2];
      STATISTICS<ComplexField_XYZ> wave[3];
      STATISTICS<ComplexField_XYZ>  pot;
      
      for (int ihad=0; ihad<2; ihad++)
         corr[ihad].input_data_bin(ibase+"/correlator."+ihad_str[ihad]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-1; dt<=+1; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            wave[dt+1].input_data_bin(ibase+"/NBSwave."+ich_str+"_"+ich_str+
                                      ".spin_1_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         pot.mem_alloc(wave[0].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0].Ndata(); i++) {
            ComplexField_XYZ tmp_Rcor[3];
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1] = (wave[dt+1](i).rot_proj(ROT_REP_A1) /
                                 (corr[0](i)(itime+dt) * corr[1](i)(itime+dt)));
            
            pot(i) = Potential::get_potential_T2(tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1]);
         }
         string ofname_base = obase+"/Pot."+ich_str+"_single.eff.3S1.t"+time_str+".bin_size"+Bsize_str;
         
         pot.output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
         pot.output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
         pot.output_data_bin       (ofname_base+".bin");
      }
      printf("END\n");
   }
   
   ////////////////////////////////////////////////////////////////
   /// For construct the leading order potential (e.g. Sigma N) ///
   {
      printf("@@@ leading order potential (for  Sigma N) calculating... "); fflush(stdout);
      
      // V_0^0, V_sigma^0, V_0^tau, V_sigma^tau
      double M4[4*4] =
      {  1, -3, -4, 12,    // V_C(I=1/2, J=0)
         1,  1, -4, -4,    // V_C(I=1/2, J=1)
         1, -3,  2, -6,    // V_C(I=3/2, J=0)
         1,  1,  2,  2  }; // V_C(I=3/2, J=1)
      
      // V_T^0, V_T^tau
      double M2[2*2] =
      {  1, -4,    // V_T(I=1/2, J=1)
         1,  2  }; // V_T(I=3/2, J=1)
      
      double M4_inv[4*4], M2_inv[2*2];
      matrix_func::inverse_matrix(M4, M4_inv, 4);
      matrix_func::inverse_matrix(M2, M2_inv, 2);
      
      string ich_str_C[4] = {
         "S"+quark+"N12__S"+quark+"N12_.V_C.J_0",
         "S"+quark+"N12__S"+quark+"N12_.V_C.J_1",
         "S"+quark+"N32__single.V_C.J_0",
         "S"+quark+"N32__single.V_C.J_1"
      };
      string ich_str_T[2] = {
         "S"+quark+"N12__S"+quark+"N12_.V_T.J_1",
         "S"+quark+"N32__single.V_T.J_1",
      };
      string och_str_C[4] = {
         "S"+quark+"N_leading.V_0_0",
         "S"+quark+"N_leading.V_s_0",
         "S"+quark+"N_leading.V_0_t",
         "S"+quark+"N_leading.V_s_t"
      };
      string och_str_T[2] = {
         "S"+quark+"N_leading.V_T_0",
         "S"+quark+"N_leading.V_T_t",
      };
      
      STATISTICS<ComplexField_XYZ> ipotC[4], ipotT[2], opotC[4], opotT[2];
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int ich=0; ich<4; ich++) {
            ipotC[ich].input_data_bin(obase+"/Pot."+ich_str_C[ich]+".t"+time_str+".bin_size"+Bsize_str+".bin");
            opotC[ich].mem_alloc(ipotC[ich].Ndata());
         }
         for (int ich=0; ich<2; ich++) {
            ipotT[ich].input_data_bin(obase+"/Pot."+ich_str_T[ich]+".t"+time_str+".bin_size"+Bsize_str+".bin");
            opotT[ich].mem_alloc(ipotT[ich].Ndata());
         }
#pragma omp parallel for
         for (int i=0; i<ipotC[0].Ndata(); i++) {
            for (int ich=0; ich<4; ich++) {
               opotC[ich](i).mem_alloc(ipotC[ich](i).get_xSIZE());
               opotC[ich](i) = 0.0;
               for (int jch=0; jch<4; jch++) opotC[ich](i) += M4_inv[4*ich+jch]*ipotC[jch](i);
            }
            for (int ich=0; ich<2; ich++) {
               opotT[ich](i).mem_alloc(ipotT[ich](i).get_xSIZE());
               opotT[ich](i) = 0.0;
               for (int jch=0; jch<2; jch++) opotT[ich](i) += M2_inv[2*ich+jch]*ipotT[jch](i);
            }
         }
         for (int ich=0; ich<4; ich++) {
            string ofname_base = obase+"/Pot."+och_str_C[ich]+".t"+time_str+".bin_size"+Bsize_str;
            
            opotC[ich].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
            opotC[ich].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
            opotC[ich].output_data_bin       (ofname_base+".bin");
         }
         for (int ich=0; ich<2; ich++) {
            string ofname_base = obase+"/Pot."+och_str_T[ich]+".t"+time_str+".bin_size"+Bsize_str;
            
            opotT[ich].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
            opotT[ich].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
            opotT[ich].output_data_bin       (ofname_base+".bin");
         }
      }
      printf("END\n");
   }
   
   /////////////////////////////////////////////////////////////////
   /// For construct the leading order potential (e.g. Lambda N) ///
   {
      printf("@@@ leading order potential (for Lambda N) calculating... "); fflush(stdout);
      
      // V_0^0, V_s^0
      double M2[2*2] =
      {  1, -3,    // V_C(I=1/2, J=0)
         1,  1  }; // V_C(I=1/2, J=1)
      
      double M2_inv[2*2];
      matrix_func::inverse_matrix(M2, M2_inv, 2);
      
      string ich_str_C[2] = {
         "L"+quark+"N12__single.V_C.J_0",
         "L"+quark+"N12__single.V_C.J_1"
      };
      string och_str_C[2] = {
         "L"+quark+"N_leading.V_0_0",
         "L"+quark+"N_leading.V_s_0",
      };
      
      STATISTICS<ComplexField_XYZ> ipotC[2], opotC[2];
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int ich=0; ich<2; ich++) {
            ipotC[ich].input_data_bin(obase+"/Pot."+ich_str_C[ich]+".t"+time_str+".bin_size"+Bsize_str+".bin");
            opotC[ich].mem_alloc(ipotC[ich].Ndata());
         }
#pragma omp parallel for
         for (int i=0; i<ipotC[0].Ndata(); i++) {
            for (int ich=0; ich<2; ich++) {
               opotC[ich](i).mem_alloc(ipotC[ich](i).get_xSIZE());
               opotC[ich](i) = 0.0;
               for (int jch=0; jch<2; jch++) opotC[ich](i) += M2_inv[2*ich+jch]*ipotC[jch](i);
            }
         }
         for (int ich=0; ich<2; ich++) {
            string ofname_base = obase+"/Pot."+och_str_C[ich]+".t"+time_str+".bin_size"+Bsize_str;
            
            opotC[ich].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
            opotC[ich].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
            opotC[ich].output_data_bin       (ofname_base+".bin");
         }
      }
      printf("END\n");
   }
   
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
