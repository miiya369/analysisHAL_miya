//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Definition of number of 1/48 reduced data points
 * @author  Takaya Miyamoto
 * @since   Thu Sep  3 01:26:45 JST 2015
 */
//--------------------------------------------------------------------------

#include <AnalysisHAL.h>

int anaHAL::reduced_Ndata(const int xSIZE, const int ySIZE, const int zSIZE) {
   int ret = 0;
   if (xSIZE == ySIZE && ySIZE == zSIZE) {
      if      (xSIZE == 16) ret = 165;
      else if (xSIZE == 32) ret = 969;
      else if (xSIZE == 48) ret = 2925;
      else if (xSIZE == 64) ret = 6545;
      else if (xSIZE == 80) ret = 12341;
      else if (xSIZE == 96) ret = 20825;
   }
   if (ret == 0) {
      for (      int z=0; z<zSIZE; z++)
         for (   int y=0; y<ySIZE; y++)
            for (int x=0; x<xSIZE; x++) {
               if (x>xSIZE/2.0 || y>x || z>y) continue;
               ret++;
            }
   }
   return ret;
}
