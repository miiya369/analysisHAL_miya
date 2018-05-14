#include <NBSwave/NBSwave.h>

static int         Lsize = 16;
static int         Tsize = 32;
static double LatSpacing = 0.1215;

static int       Dwave_m = 0;
static double factor_ang = 1.0;

static void       usage(const char*);
static void output_wave(const NBS_WAVE&, const NBS_WAVE&);

int main(int argc, char** argv)
{ 
   if (argc     == 1) { usage(argv[0]); return -1; }
   if (argc % 2 == 0) { ERROR_COMMENTS("#.arguments should be an even number."); }
   
   int n_add = (argc - 1) / 2;
   
   double* factors = new double[n_add];
   string*  ifiles = new string[n_add];
   
   for (int i=0; i<n_add; i++) {
      factors[i] = atof(argv[2*i+1]);
      ifiles [i] =      argv[2*i+2] ;
   }
   
   analysis::set_global_params(Lsize, Tsize, LatSpacing);
   
   NBS_WAVE wave, tmp_wave;
   for (int n=0; n<wave.data_size(); n++) wave(n) = COMP_ZERO;
   
   for (int i=0; i<n_add; i++) {
      tmp_wave.input_FromPath(ifiles[i].c_str());
      wave += factors[i] * tmp_wave;
   }
   
   NBS_WAVE wave_out_1(wave);
   NBS_WAVE wave_out_2(wave);
   
   // Change here ==================
   NBS_WAVE nbs_tmp1(wave);
   NBS_WAVE nbs_tmp2(wave);
   
   NBSwave::Swave_division(wave_out_1, wave_out_2);
   
   NBSwave::remove_angular(wave_out_2, Dwave_m, factor_ang);
   NBSwave::Swave_projection(wave_out_2);
   
   NBSwave::angmom_projection(nbs_tmp1, ROT_REP_E);
   NBSwave::angmom_projection(nbs_tmp2, ROT_REP_T2);
   
   wave_out_1 = nbs_tmp1 + nbs_tmp2;
   
   NBSwave::remove_angular(wave_out_1, Dwave_m, factor_ang);
   //NBSwave::angmom_projection(wave_out_1, ROT_REP_A1);
   
   //NBSwave::parity_projection(wave_out_1);
   //NBSwave::parity_projection(wave_out_2);
   // ==============================
   
   //output_wave(wave_out_1, wave_out_2);

   delete [] factors;
   delete [] ifiles;
   
   return 0;
}

static void usage(const char* argv0)
{
   printf("\nusage: %s [factor1] [ifile1] [factor2] [ifile2] ...\n\n"
	  "Then, (factor1)*(ifile1) + (factor2)*(ifile2) + ... will be calculated,\n"
	  "and output the S-wave components and those remaining (it may be D-wave components)\n"
	  "as a function of radius. (Output format: r Swave.real Swave.imag Dwave.real Dwave.imag)\n\n",
	  argv0);
}

static void output_wave(const NBS_WAVE& wave, const NBS_WAVE& wave_bar)
{
   printf("# Output format: r : Swave.real : Swave.imag : Dwave.real : Dwave.imag\n#\n");
   double R;
   for       (int iz=0; iz<Lsize/2; iz++)
      for    (int iy=0; iy<Lsize/2; iy++)
         for (int ix=0; ix<Lsize/2; ix++) {
	    
	    // If you want to print all points, comments out this line.
	    //if (ix>Lsize/2 || iy>ix || iz>iy) continue;
	    // ======================================================
	    
	    R = sqrt(ix*ix + iy*iy + iz*iz) * LatSpacing;
	    
	    printf("%lf %e %e %e %e\n", R,
		   wave    (ix,iy,iz).real(), wave    (ix,iy,iz).imag(),
		   wave_bar(ix,iy,iz).real(), wave_bar(ix,iy,iz).imag());
	 }
}
