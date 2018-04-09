//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup CORE_util_group
 * @ingroup Header_File_CORE
 * @brief   This is for implementing four point functions which consist of
 *          four Pause spinors.
 *
 * @author N.Ishii
 * @since  Thu Aug 12 23:21:16 JST 2010
 *
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup CORE_util_group CORE util
 * @ingroup  Index_group
 *
 *
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_PPPP_H
#define IS_INCLUDED_PPPP_H

#include <complex>
#include <yukawa/yukawa.h>
#include <yukawa/GPM.h>
#include <yukawa/matrix.h>
#include <yukawa/misc.h>
#include <yukawa/cubic_group.h>
#include <yukawa.CORE/PP.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

//--------------------------------------------------------------------------
/**
 * @brief Collection of functions, which may be convenient
 *        to implement PPPP
 */
//--------------------------------------------------------------------------

namespace PPPP_util {
  typedef std::complex<double> Complex;

  extern void initialize();
  extern void rotate(            Complex* buf, int ith);
  extern void charge_conjugation(Complex* buf);
  extern void cps2standard(      Complex* buf);

  extern void read( Complex* buf, FILE* fp, bool is_change_endian = false);
  extern void write(const
		    Complex* buf, FILE* fp, bool is_change_endian = false);

  extern GPM<double>     C;
  extern GPM<double>     CxC;
  extern GPM<double>     CxCxCxC;

  extern matrix<Complex> S[25];
  extern matrix<Complex> SxS[25];
  extern matrix<Complex> SstarxSstar[25];
  extern matrix<Complex> SxSxSstarxSstar[25];

  inline void copy(Complex* dst, Complex* src)
  {
    for(int i = 0; i < 2 * 2*2 * 2*2; i++)
      ((double*)dst)[i] = ((double*)src)[i];
  }
  
  inline void copy(Complex* dst, Complex* src, double factor)
  {
    for(int i = 0; i < 2 * 2*2 * 2*2; i++)
      ((double*)dst)[i] = factor * ((double*)src)[i];
  }
};

//--------------------------------------------------------------------------
/**
 * @brief sheared part of the core part for four point functions
 *        involving four Pauli spinors
 *
 *
 */
//--------------------------------------------------------------------------

template <int Length>
struct PPPP {
  typedef std::complex<double> Complex;

  Complex a[2*2 * 2*2 * Length];

  ~PPPP();
  PPPP();
  PPPP(const PPPP<Length>& other);

  PPPP<Length>& operator=(const PPPP<Length>& other);

  PPPP<Length>& operator+=(const PPPP<Length>& other);
  PPPP<Length>& operator-=(const PPPP<Length>& other);
  PPPP<Length>& unary_minus();

  PPPP<Length>& operator*=(double x);
  PPPP<Length>& operator*=(Complex z);
  PPPP<Length>& operator*=(const PPPP<Length>& other);

  PPPP<Length>& operator/=(double x);
  PPPP<Length>& operator/=(Complex z);

  const
  Complex& operator()(int alphaF,int betaF,
		      int alphaI,int betaI, int xtype) const
  { return a[range(alphaF,betaF, alphaI,betaI, xtype)]; }

  Complex& operator()(int alphaF,int betaF,
		      int alphaI,int betaI, int xtype)
  { return a[range(alphaF,betaF, alphaI,betaI, xtype)]; }

  const
  Complex& operator[](int i) const { return a[range(i)]; }
  Complex& operator[](int i)       { return a[range(i)]; }

  PPPP<Length>& rotate(int ith);
  PPPP<Length>& reflection();
  PPPP<Length>& charge_conjugation();
  PP            to_PP(Complex eta[]) const;

  void          read( FILE* fp, bool is_change_endian = false);
  void          write(FILE* fp, bool is_change_endian = false) const;
  PPPP<Length>& change_endian();

  /// This is related to N.Ishii, H.Nemura, and K.Murano
  PPPP<Length>& cps2standard();

  //--------------------------------------------------------------------------
  // static members
  //--------------------------------------------------------------------------

  /// data length in byte
  static int data_length() { return 2*2 * 2*2 * Length * sizeof(Complex); }
  static int nxtype()      { return Length; }
  static void initialize()
  {
    yukawa::cubic_group::initialize();
    PPPP_util::initialize();
  }

private:
  int range(int alphaF,int betaF, int alphaI,int betaI, int xtype) const;
  int range(int i) const;
};

/// range checkker
template <int Length>
inline int
PPPP<Length>::range(int alphaF,int betaF,int alphaI,int betaI,int xtype) const
{
  if(!NDEBUG){
    if(alphaF < 0 || 2 <= alphaF ||
       betaF  < 0 || 2 <= betaF  ||
       alphaI < 0 || 2 <= alphaI ||
       betaI  < 0 || 2 <= betaI  ||
       xtype  < 0 || Length <= xtype)
      {
	fprintf(stderr,
		"ERROR(PPPP<Length=%d>::range(int,int,int,int,int)): "
		"(%d,%d,%d,%d,%d)\n",
		Length,
		alphaF,betaF, alphaI,betaI,
		xtype);
	exit(1);
      }
  }
  return alphaF + 2*(betaF + 2*(alphaI + 2*(betaI + 2*xtype)));
}

/// range checker(2)
template <int Length>
inline int PPPP<Length>::range(int i) const
{
  if(!NDEBUG){
    if(i < 0 || 2*2 * 2*2 * Length <= i){
      fprintf(stderr,
	      "ERROR(PPPP<Length=%d>::range(int)): %d\n",
	      Length, i);
      exit(1);
    }
  }
  return i;
}

/// the destructor
template <int Length>
inline PPPP<Length>::~PPPP() {}

/// the default constructor
template <int Length>
inline PPPP<Length>::PPPP()
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] = 0.0;
}

/// the copy constructor
template <int Length>
inline PPPP<Length>::PPPP(const PPPP<Length>& other)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] = ((double*)other.a)[i];
}

/// the assignment operator
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator=(const PPPP<Length>& other)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] = ((double*)other.a)[i];
  return *this;
}


/// addition
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator+=(const PPPP<Length>& other)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] += ((double*)other.a)[i];
  return *this;
}

/// subtraction
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator-=(const PPPP<Length>& other)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] -= ((double*)other.a)[i];
  return *this;
}

/// unary minus
template <int Length>
inline PPPP<Length>& PPPP<Length>::unary_minus()
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] = -((double*)a)[i];
  return *this;
}


/// scalar multiplication (double)
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator*=(double x)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] *= x;
  return *this;
}

/// scalar multiplication (Complex)
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator*=(Complex z)
{
  for(int i = 0; i < 2*2 * 2*2 * Length; i++)
    ((Complex*)a)[i] *= z;
  return *this;
}

/// component-wise multiplication
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator*=(const PPPP<Length>& other)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] *= ((double*)other.a)[i];
  return *this;
}


/// scalar division (double)
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator/=(double x)
{
  for(int i = 0; i < 2 * 2*2 * 2*2 * Length; i++)
    ((double*)a)[i] /= x;
  return *this;
}

/// scalar division (Complex)
template <int Length>
inline PPPP<Length>& PPPP<Length>::operator/=(Complex z)
{
  for(int i = 0; i < 2*2 * 2*2 * Length; i++)
    ((Complex*)a)[i] /= z;
  return *this;
}


/// rotation
template <int Length>
inline PPPP<Length>& PPPP<Length>::rotate(int ith)
{
  for(int xtype = 0; xtype < Length; xtype++)
    PPPP_util::rotate(&a[2*2 * 2*2 * xtype], ith);
  return *this;
}

/// reflection
template <int Length>
inline PPPP<Length>& PPPP<Length>::reflection()
{
//   for(int xtype = 0; xtype < Length; xtype++)
//     PPPP_util::reflection(&a[2*2 * 2*2 * xtype]);
  return *this;
}

/// charge conjugation
template <int Length>
inline PPPP<Length>& PPPP<Length>::charge_conjugation()
{
  for(int xtype = 0; xtype < Length; xtype++)
    PPPP_util::charge_conjugation(&a[2*2 * 2*2 * xtype]);
  return *this;
}

/// read
template <int Length>
inline void PPPP<Length>::read( FILE* fp, bool is_change_endian)
{
  for(int xtype = 0; xtype < Length; xtype++)
    PPPP_util::read(&a[2*2 * 2*2 * xtype], fp, is_change_endian);
}

/// write
template <int Length>
inline void PPPP<Length>::write(FILE* fp, bool is_change_endian) const
{
  for(int xtype = 0; xtype < Length; xtype++)
    PPPP_util::write(&a[2*2 * 2*2 * xtype], fp, is_change_endian);
}

/// change endian
template <int Length>
inline PPPP<Length>& PPPP<Length>::change_endian()
{
  yukawa::change_endian((double*)a, 2 * 2*2 * 2*2 * Length);
  return *this;
}

/// cps2standard
template <int Length>
inline PPPP<Length>& PPPP<Length>::cps2standard()
{
  for(int xtype = 0; xtype < Length; xtype++)
    PPPP_util::cps2standard(&a[2*2 * 2*2 * xtype]);
  return *this;
}

/// PPPP<Length> to PP
template <int Length>
inline PP PPPP<Length>::to_PP(Complex eta[]) const
{
  PP ret;
  for(  int beta  = 0; beta  < 2; beta++)
    for(int alpha = 0; alpha < 2; alpha++){
      Complex sum = 0.0;
      for(    int xtype  = 0; xtype  < Length; xtype++)
	for(  int betaP  = 0; betaP  < 2;      betaP++)
	  for(int alphaP = 0; alphaP < 2;      alphaP++){
	    sum
	      += (*this)(alpha,beta, alphaP,betaP, xtype)
	      *  eta[alphaP + 2*(betaP + 2*xtype)];
	  }
      ret(alpha,beta) = sum;
    }
  return ret;
}

End_Name_Space_Yukawa

#endif
