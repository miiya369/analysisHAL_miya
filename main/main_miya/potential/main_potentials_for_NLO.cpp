#include <ComplexField_Sub.h>
#include <Potential.h>

#define Nch 1

string ifname_nbs(string opr_str, string ich_str, string jch_str, int itime, int ibin_size) {
   char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", ibin_size); string Bsize_str(Bsize_c);
   char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",     itime); string  time_str( time_c);
   return ("/home/miiya369/data/yku/48x48/analysis_20180521/LN/"+opr_str+
           "/bin/NBSwave."+opr_str+"."+ich_str+"_"+jch_str+".allspin.t"+time_str+".bin_size"+Bsize_str);
}
string ifname_corr(string opr_str, string had_str, int ibin_size) {
   char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", ibin_size); string Bsize_str(Bsize_c);
   return ("/home/miiya369/data/yku/48x48/analysis_20180521/LN/"+opr_str+
           "/bin/correlator."+opr_str+"."+had_str+".bin_size"+Bsize_str);
}
string ofname(string opr_str, string dir_name, string fhead, string ch_str, int itime, int ibin_size) {
   char Bsize_c[8]; snprintf(Bsize_c, sizeof(Bsize_c), "%02d", ibin_size); string Bsize_str(Bsize_c);
   char  time_c[8]; snprintf( time_c, sizeof( time_c), "%02d",     itime); string  time_str( time_c);
   return ("/home/miiya369/data/yku/48x48/analysis_20180521/LN/"+opr_str+
           "/"+dir_name+"/"+fhead+"."+opr_str+"."+ch_str+".t"+time_str+".bin_size"+Bsize_str);
}

///////// =========================== MAIN PART =========================== /////////
int main() {
   time_t stime, etime; time(&stime);
   
   /// Set parameterschar ///
   double lattice_spacing = 0.0907;
   int    bin_size        = 40;
   int    min_time        = 6;
   int    max_time        = 17;
   
   double mass_nuc = 0.601;
   double mass_lam = 0.631;
   double mass_sig = 0.642;
   
   double ZP_lam = 0.0448;
   double ZP_sig = 0.0448;
   
   string opr[2] = {"wall", "expo"};
   string spn[2] = { "1S0",  "3S1"};
   
   string ihad_str[2][2] = {
      {"proton_CG05_CG05", "proton_CG05_CG05"},
      {"Lambda_CG05_CG05", "Sigma_CG05_CG05"},
   };
   string ich_str[3] = {"LamN12", "SigN12", "SigN32"};
   double mass[2][2] = {
      {mass_nuc, mass_nuc},
      {mass_lam, mass_sig}
   };
   cdouble ZP[2] = {ZP_lam, ZP_sig};
   
   for (int iopr=0; iopr<2; iopr++)
   {
      STATISTICS<ComplexField_T> corr[2][Nch];
      for (   int ihad=0; ihad<  2; ihad++)
         for (int  ich=0;  ich<Nch; ich++ )
            corr[ihad][ich].input_data_bin(ifname_corr(opr[iopr], ihad_str[ihad][ich], bin_size));
      
      for (int ispin=0; ispin<2; ispin++)
      {
         STATISTICS<ComplexField_XYZ>  wave[3][Nch][Nch];
         STATISTICS<ComplexField_XYZ>   pot   [Nch][Nch];
         STATISTICS<ComplexField_AXYZB> wave_in;
         
         for (int itime=min_time; itime<=max_time; itime++)
         {
            printf("@@@ opr=%s, spin=%s, t=%d calculating... ", opr[iopr].c_str(), spn[ispin].c_str(), itime); fflush(stdout);
            
            for (int dt=-1; dt<=+1; dt++) {
               for (   int ich=0; ich<Nch; ich++)
                  for (int jch=0; jch<Nch; jch++) {
                     wave_in.input_data_bin(ifname_nbs(opr[iopr], ich_str[ich], ich_str[jch], itime+dt, bin_size));
                     wave[dt+1][ich][jch].mem_alloc(wave_in.Ndata());
                     for (int i=0; i<wave_in.Ndata(); i++) {
                        if (ispin == 0)
                           wave[dt+1][ich][jch](i)  = wave_in(i).spin_proj(HH_OctOct, 0, 0, HH_OctOct, 0, 0);
                        else {
                           wave[dt+1][ich][jch](i)  = wave_in(i).spin_proj(HH_OctOct, 1,+1, HH_OctOct, 1,+1);
                           wave[dt+1][ich][jch](i) += wave_in(i).spin_proj(HH_OctOct, 1, 0, HH_OctOct, 1, 0);
                           wave[dt+1][ich][jch](i) += wave_in(i).spin_proj(HH_OctOct, 1,-1, HH_OctOct, 1,-1);
                           wave[dt+1][ich][jch](i) /= 3.0;
                        }
                        wave[dt+1][ich][jch](i) = (wave[dt+1][ich][jch](i).rot_proj(ROT_REP_A1) /
                                                   (corr[0][ich](i)(itime+dt) * corr[1][ich](i)(itime+dt)));
                     }
                  }
            }
            for (   int ich=0; ich<Nch; ich++)
               for (int jch=0; jch<Nch; jch++) pot[ich][jch].mem_alloc(wave[0][ich][jch].Ndata());
            
#pragma omp parallel for
            for (int i=0; i<wave[0][0][0].Ndata(); i++) {
               ComplexField_XYZ tmp_Rcor[3][Nch*Nch], tmp_pot[Nch*Nch];
               cdouble          tmp_fac [Nch*Nch];
               
               for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
                  int ijch = jch + Nch * ich;
                  tmp_fac[ijch] = ((ZP[ich] * corr[0][jch](i)(itime) * corr[1][jch](i)(itime)) /
                                   (ZP[jch] * corr[0][ich](i)(itime) * corr[1][ich](i)(itime)));
                  
                  for (int dt=-1; dt<=+1; dt++)
                     tmp_Rcor[dt+1][ijch] = wave[dt+1][ich][jch](i);
               }
               if (Nch == 1)
                  pot[0][0](i) = Potential::get_potential_T2(tmp_Rcor[0][0], tmp_Rcor[1][0], tmp_Rcor[2][0], mass[0][0], mass[1][0]);
               else {
                  Potential::calc_CCpotential_T2(tmp_pot, tmp_Rcor[0], tmp_Rcor[1], tmp_Rcor[2], mass[0], mass[1], tmp_fac, 2);
                  
                  for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
                     int ijch = jch + Nch * ich;
                     pot[ich][jch](i) = tmp_pot[ijch];
                  }
               }
            }
            for (int ich=0; ich<Nch; ich++) for (int jch=0; jch<Nch; jch++) {
               string ofname_base;
               if (Nch == 1)
                  ofname_base = ofname(opr[iopr], "pot", "Pot", ich_str[ich]+"_single.eff."+spn[ispin], itime, bin_size);
               else
                  ofname_base = ofname(opr[iopr], "pot", "Pot", ich_str[ich]+"_"+ich_str[jch]+".eff."+spn[ispin], itime, bin_size);
               pot[ich][jch].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
               pot[ich][jch].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
               
               if (Nch == 1)
                  ofname_base = ofname(opr[iopr], "wave", "Rcr", ich_str[ich]+"_"+ich_str[jch]+".eff."+spn[ispin], itime, bin_size);
               else
                  ofname_base = ofname(opr[iopr], "wave", "Rcr", ich_str[ich]+"_single.eff."+spn[ispin], itime, bin_size);
               wave[1][ich][jch].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
               wave[1][ich][jch].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
               
               for (int i=0; i<wave[0][0][0].Ndata(); i++)
                  wave[1][ich][jch](i) = wave[1][ich][jch](i).lap();
               if (Nch == 1)
                  ofname_base = ofname(opr[iopr], "wave", "LpR", ich_str[ich]+"_"+ich_str[jch]+".eff."+spn[ispin], itime, bin_size);
               else
                  ofname_base = ofname(opr[iopr], "wave", "LpR", ich_str[ich]+"_single.eff."+spn[ispin], itime, bin_size);
               wave[1][ich][jch].output_data_err       (ofname_base+".gnu", lattice_spacing,  true);
               wave[1][ich][jch].output_data_bin_reduce(ofname_base+".fit", lattice_spacing, false);
            }
            printf("END\n");
         }
      }
   }
   
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
