#include <common/config_template.h>
#include <potential/potential.h>

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(32, 64, 0.0907);
   CHANNEL_TYPE ch("BBwave.dir.S1.00", "CG05", "CG05");
   char txyz_shift[] = "A64.000.000.000";
   
   char base[1024], fconf_list[1024], ofname[1024];
   char ibase[]    =  "/home/miiya369/data/pacs-cs/results_kekb_pacs-cs_bridge_";
   char obase[]    =  "/home/miiya369/data/pacs-cs/analysis_new";
   char dir [2][8] = {"Clover", "RHQ"};
   char odir[2][8] = {"s", "c"};
   int  bin [3]    = {57, 40, 45};
   
   double mass1[2][3] = {{0.726840, 0.642825, 0.558515}, {0.726840, 0.642825, 0.558515}};
   double mass2[2][3] = {{0.754507, 0.686167, 0.616926}, {1.233977, 1.174375, 1.118781}};
   
   int t_min = 7, t_max = 17;
   
   int spin = SPIN_1_0;
   
   NBSwave::rot_matrix_init();
   
   for (   int isc  = 0; isc  < 2; isc ++) {
      for (int iens = 0; iens < 3; iens++) {
         snprintf(base,       sizeof(base),       "%s%s/Ens%d_bin%d", ibase, dir[isc], iens+1, bin[iens]);
         snprintf(fconf_list, sizeof(fconf_list), "%s/conf.lst",      base);
         analysis::set_confs(fconf_list);
         
         R_CORRELATOR_SRC_PRJ Rcorr_sprj[3];
         R_CORRELATOR_SRC_PRJ S12_Rcorr_sprj, K_Rcorr_sprj;
         CONFIG<R_CORRELATOR> pot[2];
         pot[0].mem_alloc(analysis::Nconf); pot[1].mem_alloc(analysis::Nconf);
         
         for (int it = t_min; it <= t_max; it++) {
            for (int i=0; i<analysis::Nconf; i++) {
               for (int dt = -1; dt < 2; dt++)
                  Rcorr_sprj[dt+1].set(base, analysis::gconf(i), ch, spin, txyz_shift, it+dt);
               
               Rcorrelator::S12_psi(Rcorr_sprj[1], S12_Rcorr_sprj);
               
               potential::potential_kernels_bare(K_Rcorr_sprj, Rcorr_sprj[0], Rcorr_sprj[1], Rcorr_sprj[2],
                                                 mass1[isc][iens], mass2[isc][iens]);
               
               potential::tensor_potential(pot[0](i), pot[1](i), Rcorr_sprj[1], S12_Rcorr_sprj, K_Rcorr_sprj, spin);
            }
            snprintf(ofname, sizeof(ofname), "%s/%s_ens%d/tensor/L%sN_central_ens%d_t%02d_err",
                     obase, odir[isc], iens+1, odir[isc], iens+1, it);
            pot[0].output_data_err(ofname, analysis::lattice_spacing, true);
            snprintf(ofname, sizeof(ofname), "%s/%s_ens%d/tensor/L%sN_central_ens%d_t%02d_bin",
                     obase, odir[isc], iens+1, odir[isc], iens+1, it);
            pot[0].output_data_bin(ofname, analysis::lattice_spacing, false);
            
            snprintf(ofname, sizeof(ofname), "%s/%s_ens%d/tensor/L%sN_tensor_ens%d_t%02d_err",
                     obase, odir[isc], iens+1, odir[isc], iens+1, it);
            pot[1].output_data_err(ofname, analysis::lattice_spacing, true);
            snprintf(ofname, sizeof(ofname), "%s/%s_ens%d/tensor/L%sN_tensor_ens%d_t%02d_bin",
                     obase, odir[isc], iens+1, odir[isc], iens+1, it);
            pot[1].output_data_bin(ofname, analysis::lattice_spacing, false);
         }
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
