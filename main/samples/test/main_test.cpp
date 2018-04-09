#include <ComplexField_Sub.h>
#include <Potential.h>

namespace {
   int test1();
   int test2();
   int test3();
   int test4();
}

int main() {
   if (test3() == 0) return  0;
   else              return -1;
}

namespace {
   int test4() {
      string idir = "/home/tkymiya/gwdata/work.pacs-cs/work.run/analysis/c/bin";
      
      STATISTICS<ComplexField_XYZ  > wave0, wave1, wave_tmp0, wave_tmp1;
      STATISTICS<ComplexField_AXYZB> waveA;
      
      wave0.input_data_bin(idir+"/NBSwave.ScN32__ScN32_.spin_0_.t10.bin_size57");
      wave1.input_data_bin(idir+"/NBSwave.ScN32__ScN32_.spin_1_.t10.bin_size57");
      waveA.input_data_bin(idir+"/NBSwave.ScN32__single.allspin.t10.bin_size57");
      
      int Nconf = wave0.Ndata();
      int Ndata = wave0(0).data_size();
      
      if (Nconf != wave0.Ndata() || Nconf != wave1.Ndata() || Nconf != waveA.Ndata()) {
         printf("#.conf is differ\n");
         return 0;
      }
      wave_tmp0.mem_alloc(Nconf);
      wave_tmp1.mem_alloc(Nconf);
      
      for (int i=0; i<Nconf; i++) {
         if (Ndata     != wave0(i).data_size() || Ndata != wave1(i).data_size() ||
             Ndata *16 != waveA(i).data_size()) {
            printf("#.data is differ\n");
            return 0;
         }
         wave_tmp0(i)  = waveA(i).spin_proj(HH_OctOct, 0, 0, HH_OctOct, 0, 0);
         wave_tmp1(i)  = waveA(i).spin_proj(HH_OctOct, 1,-1, HH_OctOct, 1,-1);
         wave_tmp1(i) += waveA(i).spin_proj(HH_OctOct, 1, 0, HH_OctOct, 1, 0);
         wave_tmp1(i) += waveA(i).spin_proj(HH_OctOct, 1,+1, HH_OctOct, 1,+1);
         wave_tmp1(i) /= 3.0;
      }
      printf("#.conf=%d, #.data=%d\n", Nconf, Ndata);
      cdouble tmp;
      double  crit = 10e-16;
      int     count = 0;
      
      for (int i=0; i<Nconf; i++) for (int n=0; n<Ndata; n++) {
         tmp = 2.0 * (wave_tmp0(i)(n) - wave0(i)(n)) / (wave_tmp0(i)(n) + wave0(i)(n));
         if (abs(tmp.real()) > crit || abs(tmp.imag()) > crit) {
            printf("beyond the crit for s=0: i=%d, n=%d, abs.r=%e, abs.i=%e\n",
                   i, n, abs(tmp.real()), abs(tmp.imag()));
            count++;
         }
         tmp = 2.0 * (wave_tmp1(i)(n) - wave1(i)(n)) / (wave_tmp1(i)(n) + wave1(i)(n));
         if (abs(tmp.real()) > crit || abs(tmp.imag()) > crit) {
            printf("beyond the crit for s=1: i=%d, n=%d, abs.r=%e, abs.i=%e\n",
                   i, n, abs(tmp.real()), abs(tmp.imag()));
            count++;
         }
      }
      printf("END: #.diff=%d/%d\n", count, Nconf*Ndata);
   }
   
   int test3() {
      ComplexField_XYZ lapYlm1 = sfunc::cfield_Ylm(1,0,32);
      ComplexField_XYZ lapYlm2 = sfunc::cfield_Ylm(1,1,32);
      
      STATISTICS<ComplexField_XYZ> test(1);
      test(0) =  lapYlm2 * lapYlm1.lap();
      test.output_data_err("./test1.gnu", 0, true);
      test(0) = (lapYlm2 * lapYlm1).lap();
      test.output_data_err("./test2.gnu", 0, true);
      
      return 0;
   }
   
   int test2() {
      for (int L=0; L<=2; L++) for (int m=-L; m<=L; m++) {
         char LM_c[8]; snprintf(LM_c, sizeof(LM_c), "l%dm%+d", L,m); string LM_s(LM_c);
         ComplexField_XYZ Ylm = sfunc::cfield_Ylm(L,m,32);
         ComplexField_XYZ lap_Ylm_eigval = Ylm / Ylm.lap();
         STATISTICS<ComplexField_XYZ> test(1);
         test(0) = lap_Ylm_eigval;
         test.output_data_err("./test."+LM_s+".gnu", 0, true);
      }
      return 0;
   }
   
   int test1() {
      STATISTICS<ComplexField_XYZ> test(1);
      
      ComplexField_AXYZB test_wave(3, 32, 2);
      
      ComplexField_XYZ   Y00 = sfunc::cfield_Ylm(0, 0,32);
      ComplexField_XYZ   Y11 = sfunc::cfield_Ylm(1,-1,32);
      ComplexField_XYZ   Y10 = sfunc::cfield_Ylm(1, 0,32);
      ComplexField_XYZ   Y20 = sfunc::cfield_Ylm(2, 0,32);
      ComplexField_XYZ   Y21 = sfunc::cfield_Ylm(2,-1,32);
      ComplexField_XYZ   Y22 = sfunc::cfield_Ylm(2, 2,32);
      for (int i=0; i<32*32*32; i++) {
         test_wave(0, i, 0) = Y00(i);
         test_wave(0, i, 1) = Y10(i);
         test_wave(1, i, 0) = Y11(i);
         test_wave(1, i, 1) = Y20(i);
         test_wave(2, i, 0) = Y21(i);
         test_wave(2, i, 1) = Y22(i);
      }
      
      for (int a=0; a<3; a++) for (int b=0; b<2; b++) {
         char tmp_c[256];
         snprintf(tmp_c, sizeof(tmp_c), "test0.a%d.b%d.gnu", a, b);
         ComplexField_XYZ tmp_wave(test_wave, a, 0, b);
         test(0) = tmp_wave;
         test.output_data_err(tmp_c,0,true);
      }
      
      test_wave.output_data_miya("./test.miya");
      test_wave.output_data_bin ("./test.bin");
      
      ComplexField_AXYZB test_wave_in(3, 32, 2);
      
      test_wave_in.input_data_bin ("./test.bin");
      
      for (int a=0; a<3; a++) for (int b=0; b<2; b++) {
         char tmp_c[256];
         snprintf(tmp_c, sizeof(tmp_c), "test1.a%d.b%d.gnu", a, b);
         ComplexField_XYZ tmp_wave(test_wave_in, a, 0, b);
         test(0) = tmp_wave;
         test.output_data_err(tmp_c,0,true);
      }
      
      ComplexField_AXYZB test_wave_in2;
      
      test_wave_in2.input_data_miya("./test.miya");
      
      for (int a=0; a<3; a++) for (int b=0; b<2; b++) {
         char tmp_c[256];
         snprintf(tmp_c, sizeof(tmp_c), "test2.a%d.b%d.gnu", a, b);
         ComplexField_XYZ tmp_wave(test_wave_in2, a, 0, b);
         test(0) = tmp_wave;
         test.output_data_err(tmp_c,0,true);
      }
      return 0;
   }
}
