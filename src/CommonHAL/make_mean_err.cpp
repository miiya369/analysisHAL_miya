//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Function to make mean and error of statistical data
 * @author  Takaya Miyamoto
 * @since   Wed Apr 27 01:19:22 JST 2016
 */
//--------------------------------------------------------------------------

#include <AnalysisHAL.h>

void anaHAL::make_mean_err(const double *data, double& mean, double& err,
                           const int Ndata, const bool is_jack_knife_data) {
   //DEBUG_LOG
   double factor = double(Ndata-1);
   if (!is_jack_knife_data) factor = 1.0 / double(Ndata-1);
   
   mean           = 0.0;
   double sq_mean = 0.0;
   for (int i=0; i<Ndata; i++) {
      mean    +=     data[i];
      sq_mean += pow(data[i], 2);
   }
   mean    /= double(Ndata);
   sq_mean /= double(Ndata);
   err      = sqrt(factor * (sq_mean - pow(mean, 2)));
}

void anaHAL::make_mean_err(const cdouble *data, cdouble& mean, cdouble& err,
                           const int Ndata, const bool is_jack_knife_data) {
   //DEBUG_LOG
   double factor = double(Ndata-1);
   if (!is_jack_knife_data) factor = 1.0 / double(Ndata-1);
   
   mean            = COMP_ZERO;
   cdouble sq_mean = COMP_ZERO;
   for (int i=0; i<Ndata; i++) {
      mean    +=        data[i];
      sq_mean += cmp_sq(data[i]);
   }
   mean    /= double(Ndata);
   sq_mean /= double(Ndata);
   err      = cmp_sqrt(factor * (sq_mean - cmp_sq(mean)));
}
