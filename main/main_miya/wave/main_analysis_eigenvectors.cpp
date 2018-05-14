#include <ComplexField_Sub.h>
#include <ComplexField_Matrix.h>

namespace  { // Definitions of local variables (The initial parameters) & local functions for main part
   int      tmin = 4;
   int      tmax = 11;
   int    Nstate = 40;
   
   string ifbase = "/xc/home/takaya.miyamoto/work/su3_16x32_alpha/analysis/LN_SN";
   
   string ifcorr[2] = {
      ifbase + "/bin/correlator.proton.CG05_CG05.alpha000.bin_size10",
      ifbase + "/bin/correlator.proton.CG05_CG05.alpha000.bin_size10"
   };
   
   string ifwave(const int it) {
      char it_c[8]; snprintf(it_c, sizeof(it_c), "%02d", it); string it_str(it_c);
      return (ifbase + "/pot.bin/NBSwave.J0.8s_8s.CG05_CG05.alpha000.t" + it_str + ".bin_size10.bin");
   }
   string ifeign(const int it) {
      char it_c[8]; snprintf(it_c, sizeof(it_c), "%02d", it); string it_str(it_c);
      return (ifbase + "/eigen/Pot.J0.8s.CG05_CG05.alpha000.t" + it_str + ".bin_size10");
   }
   
   bool set_args();
   
   STATISTICS<ComplexField_T  > corr[2];
   STATISTICS<ComplexField_XYZ> wave, evec;
   ComplexField_AT eval;
   
   void input_Rcorr  (const string, const int);
   void print_eff_E_V(const int, const int);
   void print_eff_E_R(const int, const int);
   void print_eff_A  (const int, const int);
   
} // end namespace

///////// =========================== MAIN PART =========================== /////////
int main() {
   if(!set_args()) return -1;
   
   time_t stime, etime; time(&stime);
   
   for (int ihad=0; ihad<2; ihad++) corr[ihad].input_data_bin(ifcorr[ihad].c_str());
   input_Rcorr(ifwave(tmin), tmin);
   
   print_eff_E_R(tmin, tmax); printf("\n\n");
   
   char Ns_c[8]; snprintf(Ns_c, sizeof(Ns_c), "%03d", Nstate); string Ns_str(Ns_c);
   eval.mem_alloc(Nstate, wave.Ndata());
   
   print_eff_E_V(tmin, tmax); printf("\n\n");
   
   for (int it=tmin; it<=tmax; it++) {
      printf("# <t = %02d>\n", it);
      input_Rcorr(ifwave(it), it);
      eval.input_data_text_real(ifeign(it)+".Eval.N"+Ns_str+".txt");
      
      for (int istate=0; istate<Nstate; istate++) {
         char ns_c[8]; snprintf(ns_c, sizeof(ns_c), "%03d", istate); string ns_str(ns_c);
         evec.input_data_bin(ifeign(it)+".Evec.n"+ns_str+".bin");
         
         print_eff_A(it, istate);
      }
      printf("\n\n");
   }
   
   time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
   return 0;
}
///////// =========================== MAIN PART =========================== /////////

namespace  { // Inplementations of local functions for main part
   bool set_args() {
      
      return true;
   }
   ////// ------------------------------------------------------ //////
   void input_Rcorr(const string ifname, const int it) {
      wave.input_data_bin(ifname.c_str());
      for (int i=0; i<wave.Ndata(); i++)
         wave(i) = (wave(i).rot_proj(ROT_REP_A1) /
                    (corr[0](i)(it) * corr[1](i)(it)));
   }
   ////// ------------------------------------------------------ //////
   void print_eff_E_V(const int t_min, const int t_max) {
      double *E = new double[wave.Ndata()]; 
      double  mean, err;
      
      char Ns_c[8]; snprintf(Ns_c, sizeof(Ns_c), "%03d", Nstate); string Ns_str(Ns_c);
      for (int it=t_min; it<=t_max; it++) {
	 eval.input_data_text_real(ifeign(it)+".Eval.N"+Ns_str+".txt");
         
	 printf("%d", it);
	 for (int istate=0; istate<Nstate; istate++) {
	    for (int i=0; i<wave.Ndata(); i++) E[i] = eval(istate, i).real();
	    anaHAL::make_mean_err(E, mean, err, wave.Ndata(), true);
	    printf(" %1.16e %1.16e", mean, err);
	 }
	 printf("\n");
      }
      delete [] E;
   }
   ////// ------------------------------------------------------ //////
   void print_eff_E_R(const int t_min, const int t_max) {
      input_Rcorr(ifwave(t_min), t_min);
      
      cdouble *R0 = new cdouble[wave.Ndata()];
      cdouble *R1 = new cdouble[wave.Ndata()];
      cdouble  mean, err;
      
      for (int i=0; i<wave.Ndata(); i++) {
         R0[i] = 0.0;
         for (int n=0; n<wave(i).data_size(); n++) R0[i] += wave(i)(n);
      }
      
      for (int it=t_min+1; it<=t_max; it++) {
         input_Rcorr(ifwave(it), it);
         
         for (int i=0; i<wave.Ndata(); i++) {
            R1[i] = 0.0;
            for (int n=0; n<wave(i).data_size(); n++) R1[i] += wave(i)(n);
            R0[i] = -log(R1[i]/R0[i]);
         }
         anaHAL::make_mean_err(R0, mean, err, wave.Ndata(), true);
         printf("%d %1.16e %1.16e %1.16e %1.16e\n",
                it-1, mean.real(), err.real(), mean.imag(), err.imag());
         
         for (int i=0; i<wave.Ndata(); i++) R0[i] = R1[i];
      }
      delete [] R0;
      delete [] R1;
   }
   ////// ------------------------------------------------------ //////
   void print_eff_A(const int it, const int istate) {
      cdouble *A = new cdouble[wave.Ndata()];
      cdouble  mean, err;
      ComplexField_XYZ tmp_eve;
      
      for (int i=0; i<wave.Ndata(); i++) {
         A   [i] = 0.0;
         tmp_eve = evec(i).rot_proj(ROT_REP_A1);
         
         for (int n=0; n<wave(i).data_size(); n++) A[i] += tmp_eve(n) * wave(i)(n);
         
         A[i] *= exp(eval(istate,i) * double(it));
      }
      anaHAL::make_mean_err(A, mean, err, wave.Ndata(), true);
      printf("%d %1.16e %1.16e %1.16e %1.16e\n",
             istate, mean.real(), err.real(), mean.imag(), err.imag());
      
      delete [] A;
   }
}
