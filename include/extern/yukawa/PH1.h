//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup PH1_group
 * @ingroup Header_File
 * @brief   Spinor part of baryon-baryon (four point) correlator with
 *          a single A1+ source.
 * @author  N.Ishii
 * @since   Fri May 21 21:04:15 JST 2010
 *
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup PH1_group PH1
 * @ingroup  Index_group
 * @brief    Two Pauli spinors (sink and source) with A1+ orbital source.<br>
 *           Complex   data[betaP][alphaP][beta][alpha];
 *
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup PH1_Main_Code
 * @ingroup  PH1_group Main codes which involves PH1
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_PH1_H
#define IS_INCLUDED_PH1_H

#include <stdio.h>
#include <cstdlib>
#include <complex>

#include <yukawa/yukawa.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

struct PH1 {
  typedef std::complex<double> Complex;
  
  Complex a[2*2 * 2*2];

#if 0
  ~PH1();
#endif
  PH1();
  PH1(const PH1& other);

  PH1& operator=(const PH1& other);

  PH1& operator+=(const PH1& other);
  PH1& operator-=(const PH1& other);
  PH1& unary_minus();
  
  PH1& operator*=(double x);
  PH1& operator*=(Complex x);
  PH1& operator*=(const PH1& other);

  PH1& operator/=(double x);
  PH1& operator/=(Complex x);

  const
  Complex& operator()(int alphaP,int betaP,int alpha,int beta) const;
  Complex& operator()(int alphaP,int betaP,int alpha,int beta);

  const
  Complex& operator[](int i) const { return a[range(i)]; }
  Complex& operator[](int i)       { return a[range(i)]; }

  PH1& rotate(int ith);
  PH1& reflection();
  PH1& charge_conjugation();

  PH1& read( FILE* fp, bool is_change_endian = false);
  const
  PH1& write(FILE* fp, bool is_change_endian = false) const;
  PH1& change_endian();

  static void initialize();

  /// data length in bytes
  static int data_length() { return 2*2 * 2*2 * sizeof(Complex); }

  /// This is related to N.Ishii, H.Nemura, and K.Murano
  PH1& cps2standard();

private:
  int range(int alpha, int beta, int alphaP, int betaP) const;
  int range(int i) const;
};

/// range checker
inline int PH1::range(int alpha,int beta,int alphaP,int betaP) const
{
  if(!NDEBUG){
    if(alpha  < 0 || 2 <= alpha  ||
       beta   < 0 || 2 <= beta   ||
       alphaP < 0 || 2 <= alphaP ||
       betaP  < 0 || 2 <= betaP)
      {
	fprintf(stderr,
		"ERROR(PH1::range(int,int,int,int)): (%d,%d,%d,%d)\n",
		alpha, beta, alphaP, betaP);
	exit(1);
      }
  }
  return alpha + 2*(beta + 2*(alphaP + 2*betaP));
}

/// range checker 2
inline int PH1::range(int i) const
{
  if(!NDEBUG){
    if(i < 0 || 2*2 * 2*2 <= i){
      fprintf(stderr,
	      "ERROR(PH1::range(int)):  %d\n", i);
      exit(1);
    }
  }
  return i;
}

#if 0
/// the destructor
inline PH1::~PH1() {}
#endif

/// the default constructor
inline PH1::PH1()
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++) ((double*)a)[i] = 0.0;
}

/// the copy constructor
inline PH1::PH1(const PH1& other)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++)
    ((double*)a)[i] = ((double*)other.a)[i];
}

/// the assignment operator
inline PH1& PH1::operator=(const PH1& other)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++)
    ((double*)a)[i] = ((double*)other.a)[i];
  return *this;
}

/// addition
inline PH1& PH1::operator+=(const PH1& other)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++)
    ((double*)a)[i] += ((double*)other.a)[i];
  return *this;
}

/// subtraction
inline PH1& PH1::operator-=(const PH1& other)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++)
    ((double*)a)[i] -= ((double*)other.a)[i];
  return *this;
}

/// unary minus
inline PH1& PH1::unary_minus()
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++)
    ((double*)a)[i] = -((double*)a)[i];
  return *this;
}

/// scalar multiplication (double)
inline PH1& PH1::operator*=(double x)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++) ((double*)a)[i] *= x;
  return *this;
}

/// scalar multiplication (Complex)
inline PH1& PH1::operator*=(Complex x)
{
  for(int i = 0; i < 2*2 * 2*2; i++) a[i] *= x;
  return *this;
}

/// component-size multiplication (double)
inline PH1& PH1::operator*=(const PH1& other)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++)
    ((double*)a)[i] *= ((double*)other.a)[i];
  return *this;
}

/// division (double)
inline PH1& PH1::operator/=(double x)
{
  for(int i = 0; i < 2*2 * 2*2 * 2; i++) ((double*)a)[i] /= x;
  return *this;
}

/// division (Complex)
inline PH1& PH1::operator/=(Complex x)
{
  for(int i = 0; i < 2*2 * 2*2; i++) a[i] /= x;
  return *this;
}
  
///
inline const std::complex<double>&
PH1::operator()(int alpha,int beta, int alphaP,int betaP) const
{
  return a[range(alpha,beta, alphaP,betaP)];
}

///
inline std::complex<double>&
PH1::operator()(int alpha,int beta, int alphaP,int betaP)
{
  return a[range(alpha,beta, alphaP,betaP)];
}

/// reflection
inline
PH1& PH1::reflection()
{ return *this; }

//--------------------------------------------------------------------------
// helper functions
//--------------------------------------------------------------------------

inline PH1 operator+(PH1 lhs,const PH1& rhs) { return lhs += rhs; }
inline PH1 operator-(PH1 lhs,const PH1& rhs) { return lhs -= rhs; }

inline PH1 operator*(PH1 ph1,  double x) { return ph1 *= x; }
inline PH1 operator*(double x, PH1 ph1)  { return ph1 *= x; }
inline PH1 operator/(PH1 ph1,  double x) { return ph1 /= x; }

inline PH1 operator*(PH1 ph1, std::complex<double> x) { return ph1 *= x; }
inline PH1 operator*(std::complex<double> x, PH1 ph1) { return ph1 *= x; }
inline PH1 operator/(PH1 ph1, std::complex<double> x) { return ph1 /= x; }

inline PH1 operator*(PH1 lhs,const PH1& rhs) { return lhs*=rhs; }

inline PH1 operator-( PH1 obj)          { return obj.unary_minus(); }
inline PH1 rotate(    PH1 obj, int ith) { return obj.rotate(ith); }
inline PH1 reflection(PH1 obj)          { return obj.reflection(); }
inline PH1 charge_conjugation(PH1 obj)  { return obj.charge_conjugation(); } 

extern PH1 PH1_cps2standard(PH1 obj);

End_Name_Space_Yukawa

#endif
