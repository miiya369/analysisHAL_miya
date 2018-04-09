#include <common/config_template.h>
#include <potential/potential.h>

#define Nch 2

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(48, 48, 0.0907);
   CHANNEL_TYPE ch[Nch][Nch];
   ch[0][0].set("BBwave.iso.S1.00", "CG05", "CG05");
   ch[0][1].set("BBwave.iso.S1.01", "CG05", "CG05");
   ch[1][0].set("BBwave.iso.S1.02", "CG05", "CG05");
   ch[1][1].set("BBwave.iso.S1.03", "CG05", "CG05");
   char src_name  [2][8]   = {"wall", "exp"};
   char txyz_shift[2][128] = {"B00.000.000.000", "A00.006_018_A00.A00.A00"};
   char ibase[]            =  "/home/miiya369/data/yku/48x48";
   char obase[]            =  "/home/miiya369/data/yku/48x48/ana_NLO/pot2";
   
   char ofch[2][8] = {"LN", "SN"};
   char ofsr[2][8] = {"wall", "expo"};
   
   double mass1[Nch] = {0.601638, 0.601638};
   double mass2[Nch] = {0.631774, 0.642411};
   double Zfactor    = 1.0;
   
   int spin[2] = {SPIN_0_0, SPIN_1_0};
   int t_min = 6, t_max = 17;
   
   NBSwave::rot_matrix_init();
   R_CORRELATOR Rcorr[3][Nch][Nch], tmp_pot[Nch][Nch];
   CONFIG<R_CORRELATOR> pot[Nch][Nch], Rcr[Nch][Nch], LaR[Nch][Nch];
   char base[1024], fconf_list[1024], ofname[1024];
   
   for (int isrc = 0; isrc  < 2; isrc++) {
      snprintf(base,       sizeof(base),       "%s/%s.pbc/results.sloppy.bin", ibase, src_name[isrc]);
      snprintf(fconf_list, sizeof(fconf_list), "%s/conf.lst", base);
      analysis::set_confs(fconf_list);
      for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
	 pot[i][j].mem_alloc(analysis::Nconf);
	 Rcr[i][j].mem_alloc(analysis::Nconf);
	 LaR[i][j].mem_alloc(analysis::Nconf);
      }   
      for (int ispin = 0;     ispin < 2;      ispin++) {
	 for (int it = t_min; it    <= t_max; it++) {
	    for (int ic = 0; ic < analysis::Nconf; ic++) {
	       for (int dt = -1; dt < 2; dt++) 
		  for (int i = 0; i < Nch; i++) 
		     for (int j = 0; j < Nch; j++)
		        Rcorr[dt+1][i][j].set(base, analysis::gconf(ic), ch[i][j], spin[ispin], txyz_shift[isrc], it+dt);
	       
	       potential::CCP_2x2_kernels(tmp_pot, Rcorr[0], Rcorr[1], Rcorr[2], mass1, mass2, Zfactor);
	       
	       for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
		  pot[i][j](ic) = tmp_pot[i][j];
		  Rcr[i][j](ic) = Rcorr[1][i][j];
		  LaR[i][j](ic) = Rcorr[1][i][j].lap();
	       }
	    }
	    for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
	       snprintf(ofname, sizeof(ofname), "%s/spin_%d.org/%s-%s_pot_%s_t%d.gnu",
			obase, ispin, ofch[i], ofch[j], ofsr[isrc], it);
	       pot[i][j].output_data_err(ofname, analysis::lattice_spacing, true);
	       
	       snprintf(ofname, sizeof(ofname), "%s/spin_%d.org/%s-%s_pot_%s_t%d",
			obase, ispin, ofch[i], ofch[j], ofsr[isrc], it);
	       pot[i][j].output_data_bin(ofname, analysis::lattice_spacing, false);
	       
	       snprintf(ofname, sizeof(ofname), "%s/spin_%d.org/%s-%s_Rcr_%s_t%d.gnu",
			obase, ispin, ofch[i], ofch[j], ofsr[isrc], it);
	       Rcr[i][j].output_data_err(ofname, analysis::lattice_spacing, true);
	       
	       snprintf(ofname, sizeof(ofname), "%s/spin_%d.org/%s-%s_Rcr_%s_t%d",
			obase, ispin, ofch[i], ofch[j], ofsr[isrc], it);
	       Rcr[i][j].output_data_bin(ofname, analysis::lattice_spacing, false);
	       
	       snprintf(ofname, sizeof(ofname), "%s/spin_%d.org/%s-%s_LaR_%s_t%d.gnu",
			obase, ispin, ofch[i], ofch[j], ofsr[isrc], it);
	       LaR[i][j].output_data_err(ofname, analysis::lattice_spacing, true);
	       
	       snprintf(ofname, sizeof(ofname), "%s/spin_%d.org/%s-%s_LaR_%s_t%d",
			obase, ispin, ofch[i], ofch[j], ofsr[isrc], it);
	       LaR[i][j].output_data_bin(ofname, analysis::lattice_spacing, false);
	    }
	 }
      }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
