#include <ComplexField_Sub.h>
#include <Potential.h>

namespace {
  string obase = "/Users/miiya/data/tmp_wave_spherical";
  int    Lsize = 96;
  
  void output_3D(const string ofname, const ComplexField_XYZ &owave) {
    FILE *ofp = fopen(ofname.c_str(), "w");
    for (int z=0; z<Lsize; z++) for (int y=0; y<Lsize; y++) for (int x=0; x<Lsize; x++) {
	  fprintf(ofp, "%d %d %d %e %e\n", x,y,z, owave(x,y,z).real(), owave(x,y,z).imag());
	}
    fclose(ofp);
  }
  void output_r_all(const string ofname, const ComplexField_XYZ &owave) {
    FILE *ofp = fopen(ofname.c_str(), "w");
    for (int z=0; z<Lsize; z++) for (int y=0; y<Lsize; y++) for (int x=0; x<Lsize; x++) {
	  int X = x; if (x>=Lsize/2) X = Lsize-x;
	  int Y = y; if (y>=Lsize/2) Y = Lsize-y;
	  int Z = z; if (z>=Lsize/2) Z = Lsize-z;
	  fprintf(ofp, "%lf %e %e\n", sqrt(X*X+Y*Y+Z*Z), owave(x,y,z).real(), owave(x,y,z).imag());
	}
    fclose(ofp);
  }
}

int main() {
  time_t stime, etime; time(&stime);
  
  STATISTICS<ComplexField_XYZ> wave (1);
  STATISTICS<ComplexField_XYZ> wave0(1);
  STATISTICS<ComplexField_XYZ> wave1(1);
  STATISTICS<ComplexField_XYZ> wave2(1);
  STATISTICS<ComplexField_XYZ> wave3(1);
  STATISTICS<ComplexField_XYZ> wave4(1);
  wave(0).mem_alloc(Lsize);
  
  ComplexField_XYZ Y00 = sfunc::cfield_Ylm(0, 0, Lsize);
  ComplexField_XYZ Y10 = sfunc::cfield_Ylm(1, 0, Lsize);
  ComplexField_XYZ Y20 = sfunc::cfield_Ylm(2, 0, Lsize);
  ComplexField_XYZ Y30 = sfunc::cfield_Ylm(3, 0, Lsize);
  ComplexField_XYZ Y4m4= sfunc::cfield_Ylm(4,-4, Lsize);
  ComplexField_XYZ Y40 = sfunc::cfield_Ylm(4, 0, Lsize);
  ComplexField_XYZ Y4p4= sfunc::cfield_Ylm(4,+4, Lsize);
  
  int XYZ[3], XYZsize[3] = {Lsize, Lsize, Lsize};
  for (int z=0; z<Lsize; z++) for (int y=0; y<Lsize; y++) for (int x=0; x<Lsize; x++) {
	int xyz[3] = {x,y,z};
	anaHAL::convert_origin(xyz, XYZ, XYZsize, 3);
	int X=XYZ[0], Y=XYZ[1], Z=XYZ[2];
	wave(0)(x,y,z) = -exp(-(X*X+Y*Y+Z*Z)/60.0)+3.0;
      }
  wave0(0) = Y00 * wave(0);
  wave1(0) = Y10 * wave(0);
  wave2(0) = Y20 * wave(0);
  wave3(0) = Y30 * wave(0);
  wave4(0) = Y40 * wave(0);
  
  STATISTICS<ComplexField_XYZ> owave0  (1);
  STATISTICS<ComplexField_XYZ> owave4  (1);
  STATISTICS<ComplexField_XYZ> owave02 (1);
  STATISTICS<ComplexField_XYZ> owave024(1);
  STATISTICS<ComplexField_XYZ> owaveY4m4(1);
  STATISTICS<ComplexField_XYZ> owaveY40 (1);
  STATISTICS<ComplexField_XYZ> owaveY4p4(1);
  
  STATISTICS<ComplexField_XYZ> opot_0  (1);
  STATISTICS<ComplexField_XYZ> opot_02 (1);
  STATISTICS<ComplexField_XYZ> opot_024(1);
  
  owave0  (0) = wave0(0);
  owave02 (0) = wave0(0) + wave2(0);
  owave024(0) = wave0(0) + wave2(0) + 0.02*wave4(0);
  owave4  (0) =                       0.02*wave4(0);
  owaveY4m4(0) = Y4m4;
  owaveY40 (0) = Y40;
  owaveY4p4(0) = Y4p4;
  /*
  owave0  .output_data_err(obase+"/NBSwave0__.gnu", 0, true);
  owave02 .output_data_err(obase+"/NBSwave02_.gnu", 0, true);
  owave024.output_data_err(obase+"/NBSwave024.gnu", 0, true);
  owave4  .output_data_err(obase+"/NBSwave4__.gnu", 0, true);
  owaveY4m4.output_data_err(obase+"/NBSwaveY4m4.gnu", 0, true);
  owaveY40 .output_data_err(obase+"/NBSwaveY40.gnu" , 0, true);
  owaveY4p4.output_data_err(obase+"/NBSwaveY4p4.gnu", 0, true);
  
  owave0  .output_data_bin_reduce(obase+"/NBSwave0__.bin", 0, false);
  owave02 .output_data_bin_reduce(obase+"/NBSwave02_.bin", 0, false);
  owave024.output_data_bin_reduce(obase+"/NBSwave024.bin", 0, false);
  owave4  .output_data_bin_reduce(obase+"/NBSwave4__.bin", 0, false);
  owaveY4m4.output_data_bin_reduce(obase+"/NBSwaveY4m4.bin", 0, false);
  owaveY40 .output_data_bin_reduce(obase+"/NBSwaveY40.bin" , 0, false);
  owaveY4p4.output_data_bin_reduce(obase+"/NBSwaveY4p4.bin", 0, false);
  
  owave0  (0) = owave0  (0).rot_proj(ROT_REP_A1);
  owave02 (0) = owave02 (0).rot_proj(ROT_REP_A1);
  owave024(0) = owave024(0).rot_proj(ROT_REP_A1);
  owave4  (0) = owave4  (0).rot_proj(ROT_REP_A1);
  owaveY4m4(0) = owaveY4m4(0).rot_proj(ROT_REP_A1);
  owaveY40 (0) = owaveY40 (0).rot_proj(ROT_REP_A1);
  owaveY4p4(0) = owaveY4p4(0).rot_proj(ROT_REP_A1);
  
  opot_0  (0) = owave0  (0).lap() / owave0  (0);
  opot_02 (0) = owave02 (0).lap() / owave02 (0);
  opot_024(0) = owave024(0).lap() / owave024(0);
  
  owave0  .output_data_err(obase+"/NBSwave_A1_0__.gnu", 0, true);
  owave02 .output_data_err(obase+"/NBSwave_A1_02_.gnu", 0, true);
  owave024.output_data_err(obase+"/NBSwave_A1_024.gnu", 0, true);
  owave4  .output_data_err(obase+"/NBSwave_A1_4__.gnu", 0, true);
  owaveY4m4.output_data_err(obase+"/NBSwave_A1_Y4m4.gnu", 0, true);
  owaveY40 .output_data_err(obase+"/NBSwave_A1_Y40.gnu" , 0, true);
  owaveY4p4.output_data_err(obase+"/NBSwave_A1_Y4p4.gnu", 0, true);
  
  opot_0  .output_data_err(obase+"/Pot_A1_0__.gnu", 0, true);
  opot_02 .output_data_err(obase+"/Pot_A1_02_.gnu", 0, true);
  opot_024.output_data_err(obase+"/Pot_A1_024.gnu", 0, true);
  
  owave0  .output_data_bin_reduce(obase+"/NBSwave_A1_0__.bin", 0, false);
  owave02 .output_data_bin_reduce(obase+"/NBSwave_A1_02_.bin", 0, false);
  owave024.output_data_bin_reduce(obase+"/NBSwave_A1_024.bin", 0, false);
  owave4  .output_data_bin_reduce(obase+"/NBSwave_A1_4__.bin", 0, false);
  owaveY4m4.output_data_bin_reduce(obase+"/NBSwave_A1_Y4m4.bin", 0, false);
  owaveY40 .output_data_bin_reduce(obase+"/NBSwave_A1_Y40.bin ", 0, false);
  owaveY4p4.output_data_bin_reduce(obase+"/NBSwave_A1_Y4p4.bin", 0, false);
  
  opot_0  .output_data_bin_reduce(obase+"/Pot_A1_0__.bin", 0, false);
  opot_02 .output_data_bin_reduce(obase+"/Pot_A1_02_.bin", 0, false);
  opot_024.output_data_bin_reduce(obase+"/Pot_A1_024.bin", 0, false);
  */
  output_3D(obase+"/NBSwaveY4m4.3d", Y4m4);
  output_3D(obase+"/NBSwaveY40.3d" , Y40);
  output_3D(obase+"/NBSwaveY4p4.3d", Y4p4);
  output_3D(obase+"/NBSwaveY4m4_A1.3d", Y4m4.rot_proj(ROT_REP_A1));
  output_3D(obase+"/NBSwaveY40_A1.3d" , Y40 .rot_proj(ROT_REP_A1));
  output_3D(obase+"/NBSwaveY4p4_A1.3d", Y4p4.rot_proj(ROT_REP_A1));
  //output_r_all(obase+"/NBSwaveY40.ra", Y40);
  
  time(&etime); printf("#\n# JOB END : ELAPSED TIME = %d [s]\n", (int)difftime(etime,stime));
  return 0;
}
