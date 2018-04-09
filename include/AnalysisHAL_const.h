//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Header File for definition of several enum
 * @author  Takaya Miyamoto
 * @since   Wed Nov  9 07:53:25 JST 2016
 */
//--------------------------------------------------------------------------

#ifndef ANALYSIS_HAL_CONST_H
#define ANALYSIS_HAL_CONST_H

//! enum for the data type
enum {
   DTYPE_COMPLEX_FIELD_AXYZTB,
   DTYPE_COMPLEX_FIELD_AXYZB,
   DTYPE_COMPLEX_FIELD_AXYZ,
   DTYPE_COMPLEX_FIELD_XYZB,
   DTYPE_COMPLEX_FIELD_XYZ,
   DTYPE_COMPLEX_FIELD_ATB,
   DTYPE_COMPLEX_FIELD_AT,
   DTYPE_COMPLEX_FIELD_TB,
   DTYPE_COMPLEX_FIELD_T,
};

//! enum for representaions of rotation group
enum { // Do not change the number !!
   ROT_REP_A1 = 0,
   ROT_REP_A2 = 1,
   ROT_REP_E  = 2,
   ROT_REP_T1 = 3,
   ROT_REP_T2 = 4,
};

//! enum for identifier of two-hadron system
enum {
   HH_OctOct = 0,
   HH_DecDec,
   HH_DecOct,
   HH_OctMps
};

#endif
