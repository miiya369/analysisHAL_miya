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
   
   int spin[2] = {SPIN_0_0, SPIN_1_0};
   int t_min = 6, t_max = 15;
   
   NBSwave::rot_matrix_init();
   
   for (int ibin = 0; ibin < 3; ibin++) {
     snprintf(base,       sizeof(base),       "%s%d",          ibase, bin[ibin]);
     snprintf(fconf_list, sizeof(fconf_list), "%s%d/conf.lst", obase, bin[ibin]);
     analysis::set_confs(fconf_list);
     
     R_CORRELATOR Rcorr[3];
     CONFIG<R_CORRELATOR> pot(analysis::Nconf);
     
     for (  int ispin = 0;     ispin < 2;      ispin++) {
       for (int it    = t_min; it    <= t_max; it++) {
	 for (  int i=0; i<analysis::Nconf; i++) {
	   for (int dt = -1; dt < 2; dt++)
	     Rcorr[dt+1].set(base, analysis::gconf(i), ch, spin[ispin], txyz_shift, it+dt);
           
	   potential::potential_kernels(pot(i), Rcorr[0], Rcorr[1], Rcorr[2], mass1, mass2);
	 }
	 snprintf(ofname, sizeof(ofname), "%s%d/pot/LN_spin%d_t%02d_err",
		  obase, bin[ibin], ispin, it);
	 pot.output_data_err(ofname, analysis::lattice_spacing, true);
	 
	 snprintf(ofname, sizeof(ofname), "%s%d/pot/LN_spin%d_t%02d_bin",
		  obase, bin[ibin], ispin, it);
	 pot.output_data_bin(ofname, analysis::lattice_spacing, false);
       }
     }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
