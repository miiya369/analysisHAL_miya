#include <common/config_template.h>
#include <potential/potential.h>

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(48, 48, 0.0907);
   CHANNEL_TYPE ch("BBwave.dir.S1.00", "CG05", "CG05");
   char src_name  [2][8]   = {"wall", "exp"};
   char txyz_shift[2][128] = {"B00.000.000.000", "A00.006_018_A00.A00.A00"};
   char ibase[]            =  "/home/miiya369/data/yku/48x48";
   char obase[]            =  "/home/miiya369/data/yku/48x48/ana_NLO/pot";
   
   //double mass1 = 0.601638;
   //double mass2 = 0.631774;
   
   int spin[2] = {SPIN_0_0, SPIN_1_0};
   int t_min = 6, t_max = 17;
   
   NBSwave::rot_matrix_init();
   char base[1024], fconf_list[1024], ofname[1024];
   
   for (int isrc = 0; isrc  < 2; isrc++) {
      snprintf(base,       sizeof(base),       "%s/%s.pbc/results.sloppy.bin", ibase, src_name[isrc]);
      snprintf(fconf_list, sizeof(fconf_list), "%s/conf.lst", base);
      analysis::set_confs(fconf_list);
      
      R_CORRELATOR Rcorr;
      CONFIG<R_CORRELATOR> pot(analysis::Nconf);
      
      for (int ispin = 0;     ispin < 2;      ispin++) {
	 for (int it = t_min; it    <= t_max; it++) {
 	    for (int ic = 0; ic < analysis::Nconf; ic++) {
	       Rcorr.set(base, analysis::gconf(ic), ch, spin[ispin], txyz_shift[isrc], it);
	       
	       pot(ic) = Rcorr.lap() / Rcorr;
	    }
	    snprintf(ofname, sizeof(ofname), "%s/pot_%s/LN_LRpR_spin%d_t%02d_err",
		     obase, src_name[isrc], ispin, it);
	    pot.output_data_err(ofname, analysis::lattice_spacing, true);
	    
	    snprintf(ofname, sizeof(ofname), "%s/pot_%s/LN_LRpR_spin%d_t%02d_bin",
		     obase, src_name[isrc], ispin, it);
	    pot.output_data_bin(ofname, analysis::lattice_spacing, false);
	 }
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
