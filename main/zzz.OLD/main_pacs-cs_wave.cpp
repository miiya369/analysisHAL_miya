#include <common/config_template.h>
#include <NBSwave/NBSwave.h>

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(32, 64, 0.0907);
   CHANNEL_TYPE ch("BBwave.dir.S1.00", "CG05", "CG05");
   char txyz_shift[] = "A16.000.000.000";
   
   char fconf_list[1024], ofname[1024];
   char ibase[] =  "/Users/miiya/data/results/16src/c/ens1/bin057";
   char obase[] =  "/Users/miiya/data/analysis/pacs-cs/out";
   
   int spin[2] = {SPIN_0_0, SPIN_1_p1};
   int it = 7;
   
   NBSwave::rot_matrix_init();
   
   snprintf(fconf_list, sizeof(fconf_list), "%s/conf.lst", ibase);
   analysis::set_confs(fconf_list);
   
   NBS_WAVE_ORG     nbsorg;
   NBS_WAVE_SRC_PRJ nbsprj;
   CONFIG<NBS_WAVE> nbs(analysis::Nconf);
   
   for (int ispin=0; ispin<2; ispin++) for (int jspin=0; jspin<2; jspin++) {
       for (int i=0; i<analysis::Nconf; i++) {
	 nbsorg.set(ibase, analysis::gconf(i), ch, it, txyz_shift);
	 NBSwave::spin_projection(nbsorg, nbsprj, spin[ispin]);
	 NBSwave::spin_projection(nbsprj, nbs(i), spin[jspin]);
	 NBSwave::LP_projection(nbs(i));
       }
       snprintf(ofname, sizeof(ofname), "%s/LcN_spin%d-%d_err", obase, ispin, jspin);
       nbs.output_data_err(ofname, analysis::lattice_spacing, true);
     }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
