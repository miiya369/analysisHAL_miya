#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>
#include <MatrixFunc.h>
#include <Potential.h>

#define MAX_ABS_DT 2

///////// =========================== MAIN PART =========================== /////////
int main() {
   time_t stime, etime; time(&stime);
   
   /// Set parameterschar ///
   double lattice_spacing = 0.0907;
   int    bin_size        = 57;
   int    min_time        = 8;
   int    max_time        = 14;
   
   string quark = "c";
   
   string ibase = "/home/miiya369/data/pacs-cs/analysis_20180515/"+quark+"_quark/ens1/bin";
   string obase = "/home/miiya369/data/pacs-cs/analysis_20180515/"+quark+"_quark/ens1/LxN/pot_high_tdiv";
   
   double mass_nucl  = 0.726840;
   /*
   double mass_Lam   = 0.754507;
   double mass_Sig   = 0.761660;
   */
   double mass_Lam   = 1.233977;
   double mass_Sig   = 1.277982;
   
   double sq_ZP_nucl  = 0.0638732520686627;
   /*
   double sq_ZP_Lam   = 0.0681264891598947;
   double sq_ZP_Sig   = 0.0672251624030046;
   */
   double sq_ZP_Lam   = 0.1176390522262307;
   double sq_ZP_Sig   = 0.1025945201743403;
   
   
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
      STATISTICS<ComplexField_XYZ> wave[2*MAX_ABS_DT+1][2][2];
      STATISTICS<ComplexField_XYZ>  pot[2][2];
      
      for (   int ihad=0; ihad<2; ihad++)
         for (int  ich=0;  ich<2; ich++ )
            corr[ihad][ich].input_data_bin(ibase+"/correlator."+ihad_str[ihad][ich]+".bin_size"+Bsize_str);
      
      for (int itime=min_time; itime<=max_time; itime++) {
         char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
         
         for (int dt=-MAX_ABS_DT; dt<=+MAX_ABS_DT; dt++) {
            char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt);
            string time_dt_str(time_dt_c);
            
            for (   int ich=0; ich<2; ich++)
               for (int jch=0; jch<2; jch++)
                  wave[dt+MAX_ABS_DT][ich][jch].input_data_bin(ibase+"/NBSwave."+ich_str[ich]+"_"+ich_str[jch]+
                                                      ".spin_0_.t"+time_dt_str+".bin_size"+Bsize_str);
         }
         for (   int ich=0; ich<2; ich++)
            for (int jch=0; jch<2; jch++)
               pot[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
         
#pragma omp parallel for
         for (int i=0; i<wave[0][0][0].Ndata(); i++) {
            ComplexField_XYZ tmp_Rcor[2*MAX_ABS_DT+1][2][2];
            cdouble          tmp_fac [2][2];
            
            for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
               tmp_fac[ich][jch] = ((ZP[ich] * corr[0][jch](i)(itime) * corr[1][jch](i)(itime)) /
                                    (ZP[jch] * corr[0][ich](i)(itime) * corr[1][ich](i)(itime)));
               
               for (int dt=-MAX_ABS_DT; dt<=+MAX_ABS_DT; dt++)
                  tmp_Rcor[dt+MAX_ABS_DT][ich][jch] = (wave[dt+MAX_ABS_DT][ich][jch](i).rot_proj(ROT_REP_A1) /
                                              (corr[0][ich](i)(itime+dt) * corr[1][ich](i)(itime+dt)));
            }
            
            
            
            {
               ComplexField_XYZ KR[2][2];
               
               for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
                  pot[ich][jch](i).mem_alloc(tmp_Rcor[MAX_ABS_DT][ich][jch].get_xSIZE());
                  
		  double MM =  mass[0][ich] + mass[1][ich];
                  double mu = (mass[0][ich] * mass[1][ich]) / MM;
                  double dl = (mass[0][ich] - mass[1][ich]) / MM;
                  
                  /// Change Here ///
                  KR[ich][jch] = (
				  tmp_Rcor[MAX_ABS_DT][ich][jch].lap()/(2.0*mu) +
                                  
                                  (tmp_Rcor[MAX_ABS_DT-1][ich][jch] - tmp_Rcor[MAX_ABS_DT+1][ich][jch])/ 2.0 +
                                  
                                  (tmp_Rcor[MAX_ABS_DT-1][ich][jch] + tmp_Rcor[MAX_ABS_DT+1][ich][jch] -
                                   2.0*tmp_Rcor[MAX_ABS_DT][ich][jch])*(1.0+3.0*dl*dl)/(8.0*mu) +
				   
				   (tmp_Rcor[MAX_ABS_DT+2][ich][jch] - 2.0*tmp_Rcor[MAX_ABS_DT+1][ich][jch] +
				    2.0*tmp_Rcor[MAX_ABS_DT-1][ich][jch] - tmp_Rcor[MAX_ABS_DT-2][ich][jch]) *
				   (dl*dl)/(2.0*mu*MM) +
				   
				   (tmp_Rcor[MAX_ABS_DT+2][ich][jch] - 4.0*tmp_Rcor[MAX_ABS_DT+1][ich][jch] +
                                    6.0*tmp_Rcor[MAX_ABS_DT][ich][jch] - 4.0*tmp_Rcor[MAX_ABS_DT-1][ich][jch] +
				    tmp_Rcor[MAX_ABS_DT-2][ich][jch]) *
                                   (5.0*dl*dl)/(8.0*mu*MM*MM)
                                  );
                  ///////////////////
               }
               cmatrix Rmat(2);
               cmatrix Kmat(2);
               cmatrix Vmat(2);
               
               for (int n=0; n<pot[0][0](i).data_size(); n++) {
                  for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
                     Rmat(ich,jch) = tmp_Rcor[MAX_ABS_DT][ich][jch](n);
                     Kmat(ich,jch) = KR[ich][jch](n);
                  }
                  Vmat = Kmat * Rmat.inverce();
                  
                  for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++)
                     pot[ich][jch](i)(n) = Vmat(ich,jch) / tmp_fac[ich][jch];
               }
            }
            
            
            
         }
         for (int ich=0; ich<2; ich++) for (int jch=0; jch<2; jch++) {
            string ofname_base = (obase+"/Pot."+ich_str[ich]+"_"+ich_str[jch]+
                                  ".V_C.J_0.t"+time_str+".bin_size"+Bsize_str);
            
            pot[ich][jch].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
            //pot[ich][jch].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
         }
      }
      printf("END\n");
   }
   
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
