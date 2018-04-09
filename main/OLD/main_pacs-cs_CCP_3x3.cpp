#include <common/config_template.h>
#include <potential/potential.h>

#define Nch 3

static string obase = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/tmp_results_3x3/analysis/pot";
static string fbase = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/tmp_results_3x3/results_c.sprj";
static string fbase_corr = "/gwfefs/data/G17002/miyamoto/work.pacs-cs/tmp_results_3x3/results_c.split.ave";
static string shifts = "Ave.000.000.000";
static string fspin  = "iso.1p2z+1p2.1p2z+1p2/spin1_ave";
static int t_min = 6;
static int t_max = 12;

static string ch[Nch];
static string corr_name[Nch][2];

static string path_nbs(int iconf, int ich, int jch, int time) {
  char tmp_str[8]; snprintf(tmp_str, sizeof(tmp_str), "%03d+", time);
  return (fbase + "/IsoProj.NBS_S1." + ch[ich] + "_" + ch[jch] + ".dir/" + fspin + "/" + 
	  analysis::gconf(iconf) + "/NBSwave.+" + tmp_str + shifts + "." + 
	  analysis::gconf(iconf) + ".snk_NR.src_NR");
}

int main(int argc, char** argv) {
  time_t stime, etime; time(&stime);
  
  analysis::set_global_params(32, 64, 0.0907);
  analysis::set_confs("/gwfefs/data/G17002/miyamoto/work.pacs-cs/tmp_results_3x3/conf.lst.tmp");
  
  ch[0] = "Nuc__Lam_"; ch[1] = "Nuc__Sig_"; ch[2] = "Nuc__SigS";
  corr_name[0][0] = "proton_CG05_CG05"; corr_name[0][1] = "Lambda_CG05_CG05";
  corr_name[1][0] = "proton_CG05_CG05"; corr_name[1][1] = "Sigma_CG05_CG05";
  corr_name[2][0] = "proton_CG05_CG05"; corr_name[2][1] = "Sigma32";
  
  double mass[3][2] = {{0.7253592026746957, 1.2399946429390372}, 
		       {0.7253592026746957, 1.2734188696291677}, 
		       {0.7253592026746957, 1.3094691672928473}};
  double Zfactor  = 1.0;
  
  // Define the instance
  CONFIG<NBS_WAVE>      wave(analysis::Nconf);
  CONFIG<CORRELATOR>    corr   [Nch][2];
  CONFIG<R_CORRELATOR> Rcorr[3][Nch][Nch];
  CONFIG<R_CORRELATOR>   pot   [Nch][Nch];
  for (int i=0; i<Nch; i++) {
    for (int j=0; j< 2 ; j++) corr[i][j].mem_alloc(analysis::Nconf);
    for (int j=0; j<Nch; j++) {
      pot [i][j].mem_alloc(analysis::Nconf);
      for (int k=0; k<3; k++) Rcorr[k][i][j].mem_alloc(analysis::Nconf);
    }
  }
   // Read 2pt-correlator
   for (  int i=0; i<Nch; i++) 
     for (int j=0; j< 2 ; j++) {
       for (int iconf = 0; iconf < analysis::Nconf; iconf++)
	 corr[i][j](iconf).input(fbase_corr.c_str(), analysis::gconf(iconf), "PS",
				 corr_name[i][j].c_str(), 0, shifts.c_str(), true);
       
       corr[i][j].make_JK_sample();
   }
   
   for (int it = t_min; it < t_max; it++) {
     
     // Read NBS wave
     for (int dt = -1; dt < 2; dt++)
       for (  int i=0; i<Nch; i++) 
	 for (int j=0; j<Nch; j++) {
	   for (int iconf = 0; iconf < analysis::Nconf; iconf++) {
	     string base_str = path_nbs(iconf, i, j, it+dt);
	     wave(iconf).input_FromPath(base_str.c_str());
	   }
	   wave.make_JK_sample();
	   
	   for (int iconf = 0; iconf < analysis::Nconf; iconf++) 
	     Rcorr[dt+1][i][j](iconf).set(wave(iconf), corr[i][0](iconf), corr[i][1](iconf), it+dt);
	 }
     
     for (int iconf = 0; iconf < analysis::Nconf; iconf++) {
       R_CORRELATOR Rcorr_m[Nch][Nch], Rcorr_0[Nch][Nch], Rcorr_p[Nch][Nch], tmp_pot[Nch][Nch];
       for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
	   Rcorr_m[i][j] = Rcorr[0][i][j](iconf);
	   Rcorr_0[i][j] = Rcorr[1][i][j](iconf);
	   Rcorr_p[i][j] = Rcorr[2][i][j](iconf);
	 }
       potential::CCP_3x3_kernels(tmp_pot, Rcorr_m, Rcorr_0, Rcorr_p, mass, Zfactor);
       for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) pot[i][j](iconf) = tmp_pot[i][j];
       //for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) pot[i][j](iconf) = Rcorr_0[i][j].lap();
     }   
     for (int i=0; i<Nch; i++) for (int j=0; j<Nch; j++) {
	 char tmp_str[8]; snprintf(tmp_str, sizeof(tmp_str), "%03d", it);
	 string ofname = obase + "/pot." + ch[i] + "_" + ch[j] + ".t" + tmp_str + "_err";
	 pot[i][j].output_data_err(ofname.c_str(), analysis::lattice_spacing, true);
       }
     
   } // it
   
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
