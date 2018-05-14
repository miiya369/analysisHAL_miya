#include <ComplexField_Sub.h>

int main() {
  string path_anaHAL = "/Users/miiya/Dropbox/programs/analysisHAL_miya";
  string dir_idata   = path_anaHAL + "/data/test_data/NN_wall_cp-pacs_16x32_su3_10conf";
  string dir_odata   = path_anaHAL + "/data/test_results/NN_wall_cp-pacs_16x32_su3_10conf";
  string conf_lst    = dir_idata   + "/conf.lst";
  
  string had_name = "proton_CG05_CG05";
  string nbs_name = "NBS_2Boct.Prot__Neut__Prot__Neut_";
  
  int    tSize   = 32;     // Size of time    direction
  int    sSize   = 16;     // Size of spacial direction
  int    time    = 9;      // Time slice for the potential
  int    BinSize = 2;      // Bin size for jack-knife sample
  double lat_sp  = 0.1215; // Lattice spacing [fm]
  double mass    = 1.244;  // Mass of proton  [Lattice Unit]
  
  bool input_comp_data = true; // Flag for input whether 1/48 complex data
  
  // Read gauge configuration list
  anaHAL::NameList conf(conf_lst);
  
  STATISTICS<ComplexField_T>   corr(conf.Nlist());
  STATISTICS<ComplexField_XYZ> wave[3];
  ComplexField_AXYZB wave_in(4, sSize, 4);
  
  // Read single hadron correlators
  for (int ic=0; ic<corr.Ndata(); ic++) {
    corr(ic).mem_alloc(tSize);
    corr(ic).input_data_corr(dir_idata+"/correlator.PS.dir/"+conf(ic)+"/"+
			     had_name+"_correlator.+Ave.000.000.000."+conf(ic), true);
  }
  
  // Make jack-knife samples
  corr.make_JK_sample(BinSize);
  
  for (int dt=-1; dt<=1; dt++) {
    char   time_dt_c[8]; snprintf(time_dt_c, sizeof(time_dt_c), "%+04d", time+dt);
    string time_dt_str(time_dt_c);
    
    wave[dt+1].mem_alloc(conf.Nlist());
    
    // Read NBS wave functions
    for (int ic=0; ic<wave[dt+1].Ndata(); ic++) {
      if (input_comp_data) {
	wave_in.input_data_comp(dir_idata+"/"+nbs_name+".dir/"+conf(ic)+"/NBSwave."+
				time_dt_str+"+Ave.000.000.000."+conf(ic)+".oo_CG05.oo_CG05");
	
	wave[dt+1](ic) = wave_in.spin_proj(HH_OctOct, 0, 0, HH_OctOct, 0, 0);
      }
      else {
	wave[dt+1](ic).mem_alloc(sSize);
	wave[dt+1](ic).input_data_bin(dir_idata+"/Proj."+nbs_name+".dir/spin0z+0.0z+0/"+
				      conf(ic)+"/NBSwave."+time_dt_str+"+Ave.000.000.000."+
				      conf(ic)+".oo_CG05.oo_CG05");
      }
    }
    
    // Make jack-knife samples
    wave[dt+1].make_JK_sample(BinSize);
    
    // Construct S-wave R-correlator
    for (int ic=0; ic<wave[dt+1].Ndata(); ic++)
      wave[dt+1](ic) = (wave[dt+1](ic).rot_proj(ROT_REP_A1) /
			(corr(ic)(time+dt) * corr(ic)(time+dt)));
  }
    
  STATISTICS<ComplexField_XYZ> potential(wave[1].Ndata());
  
  // Calculate potential
  for (int ic=0; ic<wave[0].Ndata(); ic++) {
    potential(ic) = (wave[1](ic).lap() / mass -
		     (wave[2](ic) - wave[0](ic)) / 2.0 +
		     (wave[0](ic) + wave[2](ic) - 2.0 * wave[1](ic)) / (4.0*mass)
		     ) / wave[1](ic);
  }
  
  // Write potential with jack-knife error
  potential.output_data_err       (dir_odata+"/potential/NN_1S0_potential.err", lat_sp, true);
  potential.output_data_bin_reduce(dir_odata+"/potential/NN_1S0_potential.bin", lat_sp, false);
  
  return 0;
}
