#include <common/config_template.h>
#include <potential/potential.h>

#define Nch 2

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(32, 64, 0.0907);
   CHANNEL_TYPE ch[Nch][Nch];
   ch[0][0].set("BBwave.iso.S1.00", "CG05", "CG05");
   ch[0][1].set("BBwave.iso.S1.01", "CG05", "CG05");
   ch[1][0].set("BBwave.iso.S1.02", "CG05", "CG05");
   ch[1][1].set("BBwave.iso.S1.03", "CG05", "CG05");
   char txyz_shift[] = "A64.000.000.000";
   
   char base[1024], fconf_list[1024], ofname[1024];
   char ibase[]    =  "/home/miiya369/data/pacs-cs/results_kekb_pacs-cs_bridge_";
   char obase[]    =  "/home/miiya369/data/pacs-cs/analysis_new";
   char dir [2][8] = {"Clover", "RHQ"};
   char odir[2][8] = {"s", "c"};
   int  bin [3]    = {57, 40, 45};
   
   double mass1[2][3][Nch] =
   {
      {{0.726840, 0.726840}, {0.642825, 0.642825}, {0.558515, 0.558515}},
      {{0.726840, 0.726840}, {0.642825, 0.642825}, {0.558515, 0.558515}}
   };
   double mass2[2][3][Nch] =
   {
      {{0.754507, 0.761514}, {0.686167, 0.699378}, {0.616926, 0.641312}},
      {{1.233977, 1.277982}, {1.174375, 1.229135}, {1.118781, 1.183502}}
   };
   double Zfactor[2][3] = {{1.00846, 1.00, 1.00}, {1.15274, 1.13273, 1.14585}};
   
   int spin[2] = {SPIN_0_0, SPIN_1_0};
   int t_min = 7, t_max = 17;
   
   NBSwave::rot_matrix_init();
   
   for (   int isc  = 0; isc  < 2; isc ++) {
      for (int iens = 0; iens < 3; iens++) {
         snprintf(base,       sizeof(base),       "%s%s/Ens%d_bin%d", ibase, dir[isc], iens+1, bin[iens]);
         snprintf(fconf_list, sizeof(fconf_list), "%s/conf.lst",      base);
         analysis::set_confs(fconf_list);
         
         R_CORRELATOR Rcorr[3][Nch][Nch], tmp_pot[Nch][Nch];
         CONFIG<R_CORRELATOR> pot[Nch][Nch];
         
         for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) pot[i][j].mem_alloc(analysis::Nconf);
         
         for (   int ispin = 0;     ispin < 2;      ispin++) {
            for (int it    = t_min; it    <= t_max; it++) {
               for (int ic=0; ic<analysis::Nconf; ic++) {
                  for (int dt = -1; dt < 2; dt++) for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++)
                     Rcorr[dt+1][i][j].set(base, analysis::gconf(ic), ch[i][j], spin[ispin], txyz_shift, it+dt);
                  
                  potential::CCP_2x2_kernels(tmp_pot, Rcorr[0], Rcorr[1], Rcorr[2],
                                             mass1[isc][iens], mass2[isc][iens], Zfactor[isc][iens]);
                  
                  for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) pot[i][j](ic) = tmp_pot[i][j];
               }
               for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
                  snprintf(ofname, sizeof(ofname), "%s/%s_ens%d/pot_2x2/L%sN_spin%d_ens%d_t%02d_CCP%d%d_err",
                           obase, odir[isc], iens+1, odir[isc], ispin, iens+1, it, i, j);
                  pot[i][j].output_data_err(ofname, analysis::lattice_spacing, true);
                  
                  snprintf(ofname, sizeof(ofname), "%s/%s_ens%d/pot_2x2/L%sN_spin%d_ens%d_t%02d_CCP%d%d_bin",
                           obase, odir[isc], iens+1, odir[isc], ispin, iens+1, it, i, j);
                  pot[i][j].output_data_bin(ofname, analysis::lattice_spacing, false);
               }
            }
         }
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
