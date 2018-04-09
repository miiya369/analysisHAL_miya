//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   The spherical harmonics as field in 3-dimension space
 * @author  Takaya Miyamoto
 * @since   Thu Sep 15 00:10:27 JST 2016
 */
//--------------------------------------------------------------------------

#include <ComplexField_Sub.h>

ComplexField_XYZ sfunc::cfield_Ylm(const int L, const int M, const int Lsize) {
   DEBUG_LOG
   typedef cdouble (*Y_FUNC)(const int, const int, const int);
   Y_FUNC  Ylm;
   
   switch (L) {
      case 0:
         switch (M) {
            case  0: Ylm = Y_0_0; break;
            default: ERROR_COMMENTS("Unknown z-component of angular momentum.");
         }; break;
      case 1:
         switch (M) {
            case -1: Ylm = Y_1_m1; break;
            case  0: Ylm = Y_1_0 ; break;
            case +1: Ylm = Y_1_p1; break;
            default: ERROR_COMMENTS("Unknown z-component of angular momentum.");
         }; break;
      case 2:
         switch (M) {
            case -2: Ylm = Y_2_p2; break;
            case -1: Ylm = Y_2_p1; break;
            case  0: Ylm = Y_2_0 ; break;
            case +1: Ylm = Y_2_m1; break;
            case +2: Ylm = Y_2_m2; break;
            default: ERROR_COMMENTS("Unknown z-component of angular momentum.");
         }; break;
      default:
         ERROR_COMMENTS("Unknown angular momentum.");
   }
   
   ComplexField_XYZ ret(Lsize);
   int XYZ[3], XYZsize[3] = {Lsize, Lsize, Lsize};
   for (      int x=0; x<Lsize; x++)
      for (   int y=0; y<Lsize; y++)
         for (int z=0; z<Lsize; z++) {
            int xyz[3] = {x,y,z};
            anaHAL::convert_origin(xyz, XYZ, XYZsize, 3);
            ret(x,y,z) = (*Ylm)(XYZ[0], XYZ[1], XYZ[2]);
         }
   return ret;
}
