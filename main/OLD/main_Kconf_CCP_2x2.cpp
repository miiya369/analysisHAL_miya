#include <common/config_template.h>
#include <potential/potential.h>

#define Nch 2

int main(int argc, char** argv) {
   time_t stime, etime; time(&stime);
   
   analysis::set_global_params(96, 96, 0.0846);
   CHANNEL_TYPE ch[Nch][Nch];
   ch[0][0].set("BBwave.iso.S1.00", "CG05", "CG05");
   ch[0][1].set("BBwave.iso.S1.01", "CG05", "CG05");
   ch[1][0].set("BBwave.iso.S1.02", "CG05", "CG05");
   ch[1][1].set("BBwave.iso.S1.03", "CG05", "CG05");
   char txyz_shift[] = "C00.000.000.000";
   
   char ibase[]      = "/NAS/HDL-ZWLC/miyamoto/results.Kconf.bin.iso/bin_size_69";
   char obase[]      = "/home/miiya369/analysisHAL_miya/results/Kconf/bin69/pot2";
   char fconf_list[] = "/home/miiya369/analysisHAL_miya/results/Kconf/bin69/conf.lst";
   char ofname[1024];
   
   double mass1[Nch] = { 0.410633, 0.410633 }; // N;       t   =[14,20],                 957.788443                MeV
   double mass2[Nch] = { 0.488122, 0.523838 }; // {L, S}; {t(L)=[16,25], t(S)=[18,25]}, {1138.529041, 1221.837564} MeV
   double Zfactor    = 1.0;
   
   int spin[2] = {SPIN_0_0, SPIN_1_0};
   int t_min = 6, t_max = 15;
   
   NBSwave::rot_matrix_init();
   analysis::set_confs(fconf_list);
   
   R_CORRELATOR Rcorr[3][Nch][Nch], tmp_pot[Nch][Nch];
   CONFIG<R_CORRELATOR> pot[Nch][Nch];
   
   for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) pot[i][j].mem_alloc(analysis::Nconf);
   
   for (  int ispin = 0;     ispin < 2;      ispin++) {
     for (int it    = t_min; it    <= t_max; it++) {
       for (int ic=0; ic<analysis::Nconf; ic++) {
	 for (int dt = -1; dt < 2; dt++) 
	   for (int i=0; i<Nch; i++) 
	     for (int j=0; j<Nch; j++)
	       Rcorr[dt+1][i][j].set(ibase, analysis::gconf(ic), ch[i][j], spin[ispin], txyz_shift, it+dt);
         
	 potential::CCP_2x2_kernels(tmp_pot, Rcorr[0], Rcorr[1], Rcorr[2], mass1, mass2, Zfactor);
	 
	 for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) pot[i][j](ic) = tmp_pot[i][j];
       }
       for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
	   snprintf(ofname, sizeof(ofname), "%s/LN_spin%d_t%02d_CCP%d%d_err", obase, ispin, it, i, j);
	   pot[i][j].output_data_err(ofname, analysis::lattice_spacing, true);
	   
	   snprintf(ofname, sizeof(ofname), "%s/LN_spin%d_t%02d_CCP%d%d_bin", obase, ispin, it, i, j);
	   pot[i][j].output_data_bin(ofname, analysis::lattice_spacing, false);
	 }
     }
   }
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
