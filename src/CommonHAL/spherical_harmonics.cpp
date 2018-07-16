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
//========================================================================//
//========================================================================//

//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(1,-1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_1_m1(const int x, const int y, const int z) {
   
   return sqrt(3.0/(8.0*PI)) * cdouble(x,-y)/sqrt(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(1,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_1_0(const int x, const int y, const int z) {
   
   return sqrt(3.0/(4.0*PI)) * z/sqrt(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(1,+1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_1_p1(const int x, const int y, const int z) {
   
   return -sqrt(3.0/(8.0*PI)) * cdouble(x,y)/sqrt(x*x+y*y+z*z);
}
//========================================================================//
//========================================================================//

//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,-2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_m2(const int x, const int y, const int z) {
   
   return sqrt(15.0/(32.0*PI)) * cdouble(x,-y)*cdouble(x,-y)/double(x*x+y*y+z*z);
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
 * @brief Function for spherical harmonics Y(2,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_0(const int x, const int y, const int z) {
   
   return sqrt(5.0/(16.0*PI)) * (2.0*z*z-x*x-y*y)/double(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,+1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_p1(const int x, const int y, const int z) {
   
   return -sqrt(15.0/(8.0*PI)) * z * cdouble(x,y)/double(x*x+y*y+z*z);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(2,+2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_2_p2(const int x, const int y, const int z) {
   
   return sqrt(15.0/(32.0*PI)) * cdouble(x,y)*cdouble(x,y)/double(x*x+y*y+z*z);
}
//========================================================================//
//========================================================================//

//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,-3)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_m3(const int x, const int y, const int z) {
   
   return sqrt(35.0/(64.0*PI)) * cdouble(x,-y)*cdouble(x,-y)*cdouble(x,-y)/pow(sqrt(x*x+y*y+z*z),3);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,-2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_m2(const int x, const int y, const int z) {
   
   return sqrt(105.0/(32.0*PI)) * z * cdouble(x,-y)*cdouble(x,-y)/pow(sqrt(x*x+y*y+z*z),3);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,-1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_m1(const int x, const int y, const int z) {
   
   return sqrt(21.0/(64.0*PI)) * (4.0*z*z-x*x-y*y) * cdouble(x,-y)/pow(sqrt(x*x+y*y+z*z),3);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_0(const int x, const int y, const int z) {
   
   return sqrt(7.0/(16.0*PI)) * (2.0*z*z*z-3.0*(x*x+y*y)*z)/pow(sqrt(x*x+y*y+z*z),3);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,+1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_p1(const int x, const int y, const int z) {
   
   return -sqrt(21.0/(64.0*PI)) * (4.0*z*z-x*x-y*y) * cdouble(x,y)/pow(sqrt(x*x+y*y+z*z),3);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,+2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_p2(const int x, const int y, const int z) {
   
   return sqrt(105.0/(32.0*PI)) * z * cdouble(x,y)*cdouble(x,y)/pow(sqrt(x*x+y*y+z*z),3);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(3,+3)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_3_p3(const int x, const int y, const int z) {
   
   return -sqrt(35.0/(64.0*PI)) * cdouble(x,y)*cdouble(x,y)*cdouble(x,y)/pow(sqrt(x*x+y*y+z*z),3);
}
//========================================================================//
//========================================================================//

//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,-4)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_m4(const int x, const int y, const int z) {
   
   return sqrt(315.0/(512.0*PI)) * cdouble(x,-y)*cdouble(x,-y)*cdouble(x,-y)*cdouble(x,-y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,-3)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_m3(const int x, const int y, const int z) {
   
   return sqrt(315.0/(64.0*PI)) * z * cdouble(x,-y)*cdouble(x,-y)*cdouble(x,-y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,-2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_m2(const int x, const int y, const int z) {
   
   return sqrt(45.0/(128.0*PI)) * (6.0*z-x*x-y*y) * cdouble(x,-y)*cdouble(x,-y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,-1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_m1(const int x, const int y, const int z) {
   
   return sqrt(45.0/(64.0*PI)) * (4.0*z*z*z-3.0*(x*x+y*y)*z) * cdouble(x,-y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,0)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_0(const int x, const int y, const int z) {
   
   int x2y2 = x*x + y*y;
   return sqrt(9.0/(256.0*PI)) * (3.0*x2y2*x2y2-24.0*x2y2*z*z+8.0*z*z*z*z)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,+1)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_p1(const int x, const int y, const int z) {
   
   return -sqrt(45.0/(64.0*PI)) * (4.0*z*z*z-3.0*(x*x+y*y)*z) * cdouble(x,y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,+2)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_p2(const int x, const int y, const int z) {
   
   return sqrt(45.0/(128.0*PI)) * (6.0*z-x*x-y*y) * cdouble(x,y)*cdouble(x,y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,+3)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_p3(const int x, const int y, const int z) {
   
   return -sqrt(315.0/(64.0*PI)) * z * cdouble(x,y)*cdouble(x,y)*cdouble(x,y)/pow(double(x*x+y*y+z*z),2);
}
//--------------------------------------------------------------------------
/**
 * @brief Function for spherical harmonics Y(4,+4)
 */
//--------------------------------------------------------------------------
cdouble sfunc::Y_4_p4(const int x, const int y, const int z) {
   
   return sqrt(315.0/(512.0*PI)) * cdouble(x,y)*cdouble(x,y)*cdouble(x,y)*cdouble(x,y)/pow(double(x*x+y*y+z*z),2);
}
