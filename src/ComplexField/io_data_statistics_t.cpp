//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Function for output the complex field data with error
 * @author  Takaya Miyamoto
 * @since   Thu Sep 15 00:10:27 JST 2016
 */
//--------------------------------------------------------------------------

#include <StatisticsTemplate.h>
#include <ComplexField_Sub.h>

//--------------------------------------------------------------------------
/**
 * @brief Function to output data with error
 */
//--------------------------------------------------------------------------
template <> void STATISTICS<ComplexField_T> ::
output_data_err(const string outfile_name,
                const double lattice_spacing,
                const bool   is_jack_knife_data) {
   DEBUG_LOG
   FILE *fp;
   if ((fp = fopen(outfile_name.c_str(), "w")) == NULL)
      ERROR_FOPEN (outfile_name.c_str());
   
   double MeVfactor = 1.0;
   char   Unit[16]  = "[Lattice Unit]";
   if (lattice_spacing != 0.0) {
      MeVfactor = hbar_c / lattice_spacing;
      snprintf(Unit, sizeof(Unit), "[MeV]");
   }
   
   cdouble *tmp_mean = new cdouble[(*this)(0).data_size()];
   cdouble *tmp_err  = new cdouble[(*this)(0).data_size()];
   
   fprintf(fp, "# t, mean.real, err.real, mean.imag, err.imag %s\n", Unit);
   (*this).make_mean_err(tmp_mean, tmp_err, is_jack_knife_data);
   
   for (int n=0; n<(*this)(0).data_size(); n++) {
      fprintf(fp, "%d %1.16e %1.16e %1.16e %1.16e\n", n,
              tmp_mean[n].real() * MeVfactor, tmp_err[n].real() * MeVfactor,
              tmp_mean[n].imag() * MeVfactor, tmp_err[n].imag() * MeVfactor);
   }
   
   delete [] tmp_mean;
   delete [] tmp_err;
   fclose(fp);
}
