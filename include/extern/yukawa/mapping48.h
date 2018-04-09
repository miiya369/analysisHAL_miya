//--------------------------------------------------------------------------
/**
 * @file
 * @brief   map involved in compressing four point correlator to 1/48 
 * @ingroup Compression
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Sun May 23 02:37:51 JST 2010
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_MAPPING48_H
#define IS_INCLUDED_MAPPING48_H

#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/4pt.h>

Start_Name_Space_Yukawa

//--------------------------------------------------------------------------
/**
 * @brief  maps for compression:
 *         ixyz(general) |-> index |-> ixyz(representative)
 * @author N.ishii
 * @since  Tue May 25 01:02:06 JST 2010
 *
 */
//--------------------------------------------------------------------------

struct mapping48 {
  struct Intger;

  typedef four_point_base::BC BC;
  typedef std::complex<double> Complex;
  typedef four_point<Intger> Intger_4pt;

  int Xsites, Ysites, Zsites;
  int XYZsites;

  int      length; //< length of the representatives
  int*     rep;    //< representatives (ix,iy,iz,0)
  Intger*  map;    //< mapping to representatives (ixyz,ith,reflection,rep_idx)
  Complex* phase;  //< phase

  ~mapping48();
  mapping48();
  mapping48(int xsites,  int ysites, int zsites);
  mapping48(const mapping48& other);
  
  mapping48& operator=(const mapping48& other);

  const
  mapping48& write(FILE* fp) const;
  mapping48& read( FILE* fp);

  const
  mapping48& write(const char* fname) const;
  mapping48& read( const char* fname);

  void construct_map_and_phase(BC xbc = four_point_base::PBC,
			       BC ybc = four_point_base::PBC,
			       BC zbc = four_point_base::PBC);

  bool is_same_mapping(const mapping48& other) const;
};

//--------------------------------------------------------------------------
/**
 * @brief 
 */
//--------------------------------------------------------------------------

struct mapping48::Intger {
  typedef std::complex<double> Complex;

  int index;
  int ith;
  int parity;
  int ixyz;
  
  ~Intger(){}
  Intger()
    : index(-1), ith(-1), parity(0), ixyz(-1) {}
  Intger(int index0, int ith0, int parity0, int ixyz0)
    : index(index0), ith(ith0), parity(parity0), ixyz(ixyz0) {}

  Intger& operator+=(const Intger& other) { return *this; }
  Intger& operator-=(const Intger& other) { return *this; }
  Intger& unary_minus()                   { return *this; }

  Intger& operator*=(double)              { return *this; }
  Intger& operator*=(Complex)             { return *this; }
  Intger& operator*=(const Intger& other) { return *this; }

  Intger& operator/=(double)              { return *this; }
  Intger& operator/=(Complex)             { return *this; }

  Intger& rotate(int ith)      { return *this; }
  Intger& reflection()         { return *this; }
  Intger& charge_conjugation() { return *this; }

  static void initialize() { cubic_group::initialize(); }
};

inline mapping48::Intger
operator-(mapping48::Intger obj) { return obj.unary_minus(); }

End_Name_Space_Yukawa

#endif
