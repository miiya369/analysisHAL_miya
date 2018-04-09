//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Matrix
 * @brief   Definition for generator matrix
 * @author  Takaya Miyamoto
 * @since   Sat Nov  4 18:09:24 JST 2017
 */
//--------------------------------------------------------------------------

#include <ComplexMatrix.h>

//--------------------------------------------------------------------------
/**
 * @brief The Pauli matrix (i = 1,2,3)
 */
//--------------------------------------------------------------------------
cmatrix generator::su2(const int i) {
   cmatrix ret(2);
   switch (i) {
      case 1:
         ret(0,0) =       0.0; ret(0,1) =        1.0;
         ret(1,0) =       1.0; ret(1,1) =        0.0;
         break;
      case 2:
         ret(0,0) =       0.0; ret(0,1) = -COMP_IMAG;
         ret(1,0) = COMP_IMAG; ret(1,1) =        0.0;
         break;
      case 3:
         ret(0,0) =       1.0; ret(0,1) =        0.0;
         ret(1,0) =       0.0; ret(1,1) =       -1.0;
         break;
      default:
         ERROR_COMMENTS("Invalid index");
   }
   return ret;
}

//--------------------------------------------------------------------------
/**
 * @brief The Gell-Mann matrix (i = 1,2,3,4,5,6,7,8)
 */
//--------------------------------------------------------------------------
cmatrix generator::su3(const int i) {
   cmatrix ret(3);
   switch (i) {
      case 1:
         ret(0,0) =       0.0; ret(0,1) =        1.0; ret(0,2) =        0.0;
         ret(1,0) =       1.0; ret(1,1) =        0.0; ret(1,2) =        0.0;
         ret(2,0) =       0.0; ret(2,1) =        0.0; ret(2,2) =        0.0;
         break;
      case 2:
         ret(0,0) =       0.0; ret(0,1) = -COMP_IMAG; ret(0,2) =        0.0;
         ret(1,0) = COMP_IMAG; ret(1,1) =        0.0; ret(1,2) =        0.0;
         ret(2,0) =       0.0; ret(2,1) =        0.0; ret(2,2) =        0.0;
         break;
      case 3:
         ret(0,0) =       1.0; ret(0,1) =        0.0; ret(0,2) =        0.0;
         ret(1,0) =       0.0; ret(1,1) =       -1.0; ret(1,2) =        0.0;
         ret(2,0) =       0.0; ret(2,1) =        0.0; ret(2,2) =        0.0;
         break;
      case 4:
         ret(0,0) =       0.0; ret(0,1) =        0.0; ret(0,2) =        1.0;
         ret(1,0) =       0.0; ret(1,1) =        0.0; ret(1,2) =        0.0;
         ret(2,0) =       1.0; ret(2,1) =        0.0; ret(2,2) =        0.0;
         break;
      case 5:
         ret(0,0) =       0.0; ret(0,1) =        0.0; ret(0,2) = -COMP_IMAG;
         ret(1,0) =       0.0; ret(1,1) =        0.0; ret(1,2) =        0.0;
         ret(2,0) = COMP_IMAG; ret(2,1) =        0.0; ret(2,2) =        0.0;
         break;
      case 6:
         ret(0,0) =       0.0; ret(0,1) =        0.0; ret(0,2) =        0.0;
         ret(1,0) =       0.0; ret(1,1) =        0.0; ret(1,2) =        1.0;
         ret(2,0) =       0.0; ret(2,1) =        1.0; ret(2,2) =        0.0;
         break;
      case 7:
         ret(0,0) =       0.0; ret(0,1) =        0.0; ret(0,2) =        0.0;
         ret(1,0) =       0.0; ret(1,1) =        0.0; ret(1,2) = -COMP_IMAG;
         ret(2,0) =       0.0; ret(2,1) =  COMP_IMAG; ret(2,2) =        0.0;
         break;
      case 8:
         ret(0,0) =       1.0; ret(0,1) =        0.0; ret(0,2) =        0.0;
         ret(1,0) =       0.0; ret(1,1) =        1.0; ret(1,2) =        0.0;
         ret(2,0) =       0.0; ret(2,1) =        0.0; ret(2,2) =       -2.0;
         ret /= sqrt(3.0);
         break;
      default:
         ERROR_COMMENTS("Invalid index");
   }
   return ret;
}
