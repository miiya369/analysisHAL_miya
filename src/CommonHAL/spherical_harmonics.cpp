//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Special functions
 * @brief   The function for spherical harmonics
 * @author  Takaya Miyamoto
 * @since   Wed Dec 16 06:20:44 JST 2015
 */
//--------------------------------------------------------------------------

#include <AnalysisHAL.h>

//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(0,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_0_0(const int x, const int y, const int z) {
   
   return 1.0/sqrt(4.0*PI);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(1,+1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_1_p1(const int x, const int y, const int z) {
   
   return -sqrt(3.0/(8.0*PI)) * cdouble(x, y)/sqrt(x*x + y*y + z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(1,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_1_0(const int x, const int y, const int z) {
   
   return sqrt(3.0/(4.0*PI)) * z/sqrt(x*x + y*y + z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(1,-1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_1_m1(const int x, const int y, const int z) {
   
   return sqrt(3.0/(8.0*PI)) * cdouble(x,-y)/sqrt(x*x + y*y + z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,+2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_p2(const int x, const int y, const int z) {
   
   return sqrt(15.0/(32.0*PI)) * cdouble(x, y)*cdouble(x, y)/double(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,+1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_p1(const int x, const int y, const int z) {
   
   return -sqrt(15.0/(8.0*PI)) * z * cdouble(x, y)/double(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_0(const int x, const int y, const int z) {
   
   return sqrt(5.0/(16.0*PI)) * (2.0*z*z - x*x - y*y)/double(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,-1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_m1(const int x, const int y, const int z) {
   
   return sqrt(15.0/(8.0*PI)) * z * cdouble(x,-y)/double(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,-2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_m2(const int x, const int y, const int z) {
   
   return sqrt(15.0/(32.0*PI)) * cdouble(x,-y)*cdouble(x,-y)/double(x*x+y*y+z*z);
}
