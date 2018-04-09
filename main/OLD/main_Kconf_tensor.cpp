#include <common/config_template.h>
#include <potential/potential.h>

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(96, 96, 0.0846);
   CHANNEL_TYPE ch("BBwave.dir.S1.00", "CG05", "CG05");
   char txyz_shift[] = "C00.000.000.000";
   
   char base[1024], fconf_list[1024], ofname[1024];
   char ibase[]    =  "/NAS/HDL-ZWLC/results.Kconf.bin/bin_size_";
   char obase[]    =  "/home/miiya369/analysisHAL_miya/results/Kconf/bin";
   int  bin [3]    = {23, 46, 69};
   
   double mass1 = 0.410633; // t=[14,20], 957.788443  MeV
   double mass2 = 0.488122; // t=[16,25], 1138.529041 MeV
   
   int t_min = 6, t_max = 15;
   
   int spin = SPIN_1_0;
   
   NBSwave::rot_matrix_init();
   
   for (int ibin = 0; ibin < 3; ibin++) {
     snprintf(base,       sizeof(base),       "%s%d",          ibase, bin[ibin]);
     snprintf(fconf_list, sizeof(fconf_list), "%s%d/conf.lst", obase, bin[ibin]);
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
         
	 potential::potential_kernels_bare(K_Rcorr_sprj, Rcorr_sprj[0], Rcorr_sprj[1], Rcorr_sprj[2], mass1, mass2);
	 
	 potential::tensor_potential(pot[0](i), pot[1](i), Rcorr_sprj[1], S12_Rcorr_sprj, K_Rcorr_sprj, spin);
       }
       snprintf(ofname, sizeof(ofname), "%s%d/tensor/LN_central_t%02d_err", obase, bin[ibin], it);
       pot[0].output_data_err(ofname, analysis::lattice_spacing, true);
       snprintf(ofname, sizeof(ofname), "%s%d/tensor/LN_central_t%02d_bin", obase, bin[ibin], it);
       pot[0].output_data_bin(ofname, analysis::lattice_spacing, false);
       
       snprintf(ofname, sizeof(ofname), "%s%d/tensor/LN_tensor_t%02d_err", obase, bin[ibin], it);
       pot[1].output_data_err(ofname, analysis::lattice_spacing, true);
       snprintf(ofname, sizeof(ofname), "%s%d/tensor/LN_tensor_t%02d_bin", obase, bin[ibin], it);
       pot[1].output_data_bin(ofname, analysis::lattice_spacing, false);
     }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
