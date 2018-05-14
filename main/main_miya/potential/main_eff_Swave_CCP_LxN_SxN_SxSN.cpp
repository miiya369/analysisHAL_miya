#include <ComplexField_Sub.h>
#include <Potential.h>
#define Nch 3

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   double lat_sp    = 0.0907;
   int bin_size     = 57;
   int min_time_pot = 10;
   int max_time_pot = 12;
   
   string quark    = "s";
   bool   is_2x2ch = false;
   
   string wdir  = "/Users/miiya/tmp/pacs-cs.fx100.LcN_3x3.bin";
   string ibase = wdir;
   string obase = wdir+"/analysis";
   
   string ihad_str[2][3] = {
      {"proton_CG05_CG05", "proton_CG05_CG05", "proton_CG05_CG05"},
      {"Lambda_CG05_CG05", "Sigma_CG05_CG05" , "Sigma32"}
   };
   double masses_s[2][3] = {
      {0.7211, 0.7211, 0.7211},
      {0.7502, 0.7568, 0.8617}
   };
   double masses_c[2][3] = {
      {0.7211, 0.7211, 0.7211},
      {1.2337, 1.2775, 1.3179}
   };
   double sq_ZP_s[2][3] = {
      {0.060580, 0.060580, 0.060580},
      {0.065147, 0.064402, 0.076991}
   };
   double sq_ZP_c[2][3] = {
      {0.060580, 0.060580, 0.060580},
      {0.117082, 0.101861, 0.141115}
   };
   
   string ch_3x3_str[3] = {"L"+quark+"N12_", "S"+quark+"N12_", "S"+quark+"SN12"};
   string ch_2x2_str[2] = {"S"+quark+"N32_", "S"+quark+"SN32"};
   
   int pad;
   STATISTICS<ComplexField_T  > corr[2][3];
   STATISTICS<ComplexField_XYZ> wave[3][Nch][Nch];
   STATISTICS<ComplexField_XYZ> pot [Nch][Nch];
   
   void calc_potential(const int);
}
/// ================================================================================= ///
int main() {
   time_t stime, etime; time(&stime);
   char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", bin_size); string Bsize_str(Bsize_c);
   string *ch_str;
   
   if (is_2x2ch) {
      ch_str = ch_2x2_str;
      pad    = 1;
   } else {
      ch_str = ch_3x3_str;
      pad    = 0;
   }
   // *** Read single hadron correlators *** //
   printf("Reading Correlators... "); fflush(stdout);
   for (int ihad=0; ihad<2; ihad++) for (int ich=0; ich<3; ich++)
      corr[ihad][ich].input_data_bin(ibase+"/corr_"+quark+"/correlator."+ihad_str[ihad][ich]+".bin_size"+Bsize_str);
   printf("end\n");
   // ************************************** //
   
   /*{
      for (int ihad=0; ihad<2; ihad++) for (int ich=0; ich<3; ich++) {
         STATISTICS<ComplexField_T> effmass(corr[ihad][ich].Ndata());
         for (int i=0; i<corr[ihad][ich].Ndata(); i++) {
            effmass(i).mem_alloc(63);
            for (int it=0; it<63; it++) effmass(i)(it) = -log(corr[ihad][ich](i)(it+1)/corr[ihad][ich](i)(it));
         } effmass.output_data_err(obase+"/Ems."+ihad_str[ihad][ich]+".bin_size"+Bsize_str+".gnu", lat_sp, true);
      }
   }*/
   for (int itime=min_time_pot; itime<=max_time_pot; itime++) {
      char time_c[8]; snprintf(time_c, sizeof(time_c), "%02d", itime); string time_str(time_c);
      
      // ********* Read 4pt-correlators ********* //
      printf("Reading NBSwave... "); fflush(stdout);
      for (int dt=-1; dt<=+1; dt++) {
         char time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%02d", itime+dt); string time_dt_str(time_dt_c);
         for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++)
            wave[dt+1][ich][jch].input_data_bin(ibase+"/NBSwave."+ch_str[ich]+"_"+ch_str[jch]+
                                                ".spin_1_.t"+time_dt_str+".bin_size"+Bsize_str);
      } printf("end\n");
      // **************************************** //
      
      calc_potential(itime);
      
      /*{
         int ich = 0, jch = ich; double m[2];
         for (int j=0; j<2; j++)
            if   (quark == "s") m[j] = masses_s[j][ich+pad];
            else                m[j] = masses_c[j][ich+pad];
         double mu = (m[0]*m[1])/(m[0]+m[1]), dl = (m[0]-m[1])/(m[0]+m[1]);
         
         STATISTICS<ComplexField_XYZ> spot(wave[1][ich][jch].Ndata());
         for (int i=0; i<wave[1][ich][jch].Ndata(); i++){
            ComplexField_XYZ R[3];
            for (int dt=-1; dt<=+1; dt++)
               R[dt+1] = (wave[dt+1][ich][jch](i).rot_proj(ROT_REP_A1) /
                          (corr[0][ich+pad](i)(itime+dt)*corr[1][ich+pad](i)(itime+dt)));
            
            spot(i) = ( R[1].lap()/(2.0*mu) +
                       (R[0] - R[2])/ 2.0 +
                       (R[0] + R[2] - 2.0*R[1])*(1.0+3.0*dl*dl)/(8.0*mu)
                       ) / R[1];
         }
         string ofbase = obase+"/Pot."+ch_str[ich]+"_single.Veff1S0.t"+time_str+".bin_size"+Bsize_str;
         spot.output_data_err       (ofbase+".gnu", lat_sp,  true);
         spot.output_data_bin_reduce(ofbase+".bin", lat_sp, false);
      }*/
      
      // ********* Write potentials ********* //
      printf("Output potentials... "); fflush(stdout);
      /*{
         int ich = 0, jch = 1;
         for (int i=0; i<wave[1][ich][jch].Ndata(); i++) {
            wave[1][ich][jch](i) =  wave[1][ich][jch](i) - wave[1][ich][jch](i).rot_proj(ROT_REP_A1);
            wave[1][ich][jch](i) = (wave[1][ich][jch](i) / (corr[0][ich+pad](i)(itime)*corr[1][ich+pad](i)(itime)));
         }
         wave[1][ich][jch].output_data_err(obase+"/Rcr."+ch_str[ich]+"_"+ch_str[jch]+".5D0.t"+
                                           time_str+".bin_size"+Bsize_str+".gnu", lat_sp, true);
      }*/
      
      for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++)
      {
         /*
         for (int i=0; i<wave[1][ich][jch].Ndata(); i++) {
            wave[1][ich][jch](i) = wave[1][ich][jch](i).rot_proj(ROT_REP_A1);
            //wave[1][ich][jch](i) = wave[1][ich][jch](i) - wave[1][ich][jch](i).rot_proj(ROT_REP_A1);
            wave[1][ich][jch](i) = (wave[1][ich][jch](i) / (corr[0][ich+pad](i)(itime)*corr[1][ich+pad](i)(itime)));
         }
         wave[1][ich][jch].output_data_err(obase+"/Rcr."+ch_str[ich]+"_"+ch_str[jch]+".3S1.t"+
                                           time_str+".bin_size"+Bsize_str+".gnu", lat_sp, true);
         */
         
         string ofbase = obase+"/Pot."+ch_str[ich]+"_"+ch_str[jch]+".Veff3S1.t"+time_str+".bin_size"+Bsize_str;
         pot[ich][jch].output_data_err       (ofbase+".gnu", lat_sp,  true);
         pot[ich][jch].output_data_bin_reduce(ofbase+".bin", lat_sp, false);
         
      } printf("end\n");
      // ************************************ //
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}

/// ================================================================================= ///

namespace {
   void calc_potential(const int itime)
   {
      double masses[2][Nch], ZP[Nch];
      for (int i=0; i<Nch; i++) {
         if (quark == "s") {
            for (int j=0; j<2; j++) masses[j][i] = masses_s[j][i+pad];
            ZP[i] = sq_ZP_s[0][i+pad] * sq_ZP_s[1][i+pad];
         } else {
            for (int j=0; j<2; j++) masses[j][i] = masses_c[j][i+pad];
            ZP[i] = sq_ZP_c[0][i+pad] * sq_ZP_c[1][i+pad];
         }
      }
      for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) pot[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
      
      printf("Calculating potentials... \n");
#pragma omp parallel for
      for (int i=0; i<wave[0][0][0].Ndata(); i++){
         ComplexField_XYZ tmp_pot[Nch*Nch], tmp_Rcor[3][Nch*Nch];
         cdouble          tmp_fac[Nch*Nch];
         
         for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
            int  ijch = jch + Nch*ich;
            tmp_fac[ijch] = ((ZP[ich] * corr[0][jch+pad](i)(itime) * corr[1][jch+pad](i)(itime)) /
                             (ZP[jch] * corr[0][ich+pad](i)(itime) * corr[1][ich+pad](i)(itime)));
            
            for (int dt=-1; dt<=+1; dt++)
               tmp_Rcor[dt+1][ijch] = (wave[dt+1][ich][jch](i).rot_proj(ROT_REP_A1) /
                                       (corr[0][ich+pad](i)(itime+dt) * corr[1][ich+pad](i)(itime+dt)));
         }
         Potential::calc_CCpotential_T2(tmp_pot, tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2],
                                        masses[0], masses[1], tmp_fac, Nch);
         
         for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
            int ijch = jch + Nch*ich;
            pot[ich][jch](i) = tmp_pot[ijch];
         } printf("Calculating potentials... end: conf=%d\n", i);
      }
   }
}
