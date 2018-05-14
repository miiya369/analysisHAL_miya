#include <ComplexField_Sub.h>

int main() {
  string path_anaHAL = "/Users/miiya/Dropbox/programs/analysisHAL_miya";
  string dir_idata   = path_anaHAL + "/data/test_data/NN_wall_cp-pacs_16x32_su3_10conf";
  string dir_odata   = path_anaHAL + "/data/test_results/NN_wall_cp-pacs_16x32_su3_10conf";
  string conf_lst    = dir_idata   + "/conf.lst";
  
  string had_name = "proton_CG05_CG05";
  
  int    tSize   = 32;     // Size of time direction
  int    BinSize = 2;      // Bin size for jack-knife sample
  double lat_sp  = 0.1215; // Lattice spacing [fm]
  
  // Read gauge configuration list
  anaHAL::NameList conf(conf_lst);
  
  STATISTICS<ComplexField_T> corr(conf.Nlist());
  
  // Read single hadron correlators
  for (int ic=0; ic<corr.Ndata(); ic++) {
    corr(ic).mem_alloc(tSize);
    corr(ic).input_data_corr(dir_idata+"/correlator.PS.dir/"+conf(ic)+"/"+
			     had_name+"_correlator.+Ave.000.000.000."+conf(ic), true);
  }
  
  // Make jack-knife samples
  corr.make_JK_sample(BinSize);
  
  STATISTICS<ComplexField_T> effmass(corr.Ndata());
  
  // Calculate effective mass
  for (int ic=0; ic<corr.Ndata(); ic++) {
    effmass(ic).mem_alloc(tSize-1);  
    for (int it=0; it<tSize-1; it++) effmass(ic)(it) = -log(corr(ic)(it+1) / corr(ic)(it));
  }
  
  // Write effective mass with jack-knife error
  effmass.output_data_err(dir_odata+"/eff_mass/"+had_name+"_effective_mass.err", lat_sp, true);
  
  return 0;
}
