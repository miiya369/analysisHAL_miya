//--------------------------------------------------------------------------
/**
 * @file
 * @brief   Core part of BS wave function of two Pauli spinors
 * @ingroup PP_group
 * @ingroup Header_File_CORE
 * @author  N.Ishii
 * @since   Sun Aug 15 18:06:25 JST 2010
 *
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup PP_group PP
 * @ingroup  Index_group
 *
 * The format:<br>
 * Complex ppwave[iz][iy][ix][beta][alpha];
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup PP_Main_Code Main codes which involves PP
 * @ingroup  PP_group
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_PP_H
#define IS_INCLUDED_PP_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/matrix.h>
#include <yukawa/GPM.h>
#include <yukawa/misc.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

struct PP {
  typedef std::complex<double> Complex;
  typedef matrix<Complex>      Matrix;

  Complex a[2*2];

  ~PP();
  PP();
  PP(const PP& other);

  PP& operator=(const PP& other);

  PP& operator+=(const PP& other);
  PP& operator-=(const PP& other);
  PP& unary_minus();

  PP& operator*=(double x);
  PP& operator*=(Complex z);
  PP& operator*=(const PP& other);

  PP& operator/=(double x);
  PP& operator/=(Complex z);

  const
  Complex& operator()(int alpha,int beta) const
  { return a[range(alpha,beta)]; }

  Complex& operator()(int alpha,int beta)
  { return a[range(alpha,beta)]; }

  const
  Complex& operator[](int i) const { return a[range(i)]; }
  Complex& operator[](int i)       { return a[range(i)]; }

  PP& rotate(int ith);
  PP& reflection();
  PP& charge_conjugation();

  void read( FILE* fp, bool is_change_endian = false);
  void write(FILE* fp, bool is_change_endian = false) const;

  PP& change_endian();

  PP& spin_singlet();
  PP& spin_triplet();

  /// data length in bytes
  static int  data_length() { return 4 * sizeof(Complex); }
  static void initialize();

  static Matrix      S[24];
  static Matrix      SxS[24];
  static GPM<double> C;
  static GPM<double> CxC;

  // additional
  enum   Rep_O {A1=0, A2=1, E=2, T2=3, T1=4, UNKNOWN=-1};
  static double character_O_ary[];
  static const
  double& character_O(int i,int rep)
  { return character_O_ary[i + 5*rep]; }

private:
  int range(int alpha,int beta) const {
    if(!NDEBUG){
      if(alpha < 0 || 2 <= alpha ||
	 beta  < 0 || 2 <= beta)
	{
	  fprintf(stderr,
		  "ERROR(PP::range(int,int)): index out of range\n");
	  exit(1);
	}
    }
    return alpha + 2*beta;
  }
  int range(int i) const {
    if(!NDEBUG){
      if (i < 0 || 4 <= i){
	fprintf(stderr,
		"ERROR(PP::range(int)): index out of range\n");
	exit(1);
      }
    }
    return i;
  }
};

/// The destructor
inline PP::~PP()
{}

/// The default constructor
inline PP::PP()
{ for(int i = 0; i < 2 * 4; i++) ((double*)a)[i] = 0.0; }

/// The copy constructor
inline PP::PP(const PP& other)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] = ((double*)other.a)[i];
}

/// The assignment operator
inline PP& PP::operator=(const PP& other)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] = ((double*)other.a)[i];
  return *this;
}

/// addition
inline PP& PP::operator+=(const PP& other)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] += ((double*)other.a)[i];
  return *this;
}

/// subtraction
inline PP& PP::operator-=(const PP& other)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] -= ((double*)other.a)[i];
  return *this;
}

/// unary minus
inline PP& PP::unary_minus()
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] = -((double*)a)[i];
  return *this;
}

/// scalar multiplication (double)
inline PP& PP::operator*=(double x)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] *= x;
  return *this;
}

/// scalar multiplication (Complex)
inline PP& PP::operator*=(Complex z)
{
  for(int i = 0; i < 4; i++)
    ((Complex*)a)[i] *= z;
  return *this;
}

/// component-wise multiplication
inline PP& PP::operator*=(const PP& other)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] *= ((double*)other.a)[i];
  return *this;
}


/// scalar division (double)
inline PP& PP::operator/=(double x)
{
  for(int i = 0; i < 2 * 4; i++)
    ((double*)a)[i] /= x;
  return *this;
}

/// scalar division (Complex)
inline PP& PP::operator/=(Complex z)
{
  for(int i = 0; i < 4; i++)
    ((Complex*)a)[i] /= z;
  return *this;
}

/// reflection
inline PP& PP::reflection()
{
  return *this;
}

/// read
inline void PP::read(FILE* fp, bool is_change_endian)
{
  int length = fread(a, sizeof(Complex), 4, fp);
  if(length != 4){
    fprintf(stderr, "ERROR(PP::read): length = %d != 4\n", length);
    exit(1);
  }
  if(is_change_endian == true){
    yukawa::change_endian((double*)a, 2 * 4);
  }
}

/// write
inline void PP::write(FILE* fp, bool is_change_endian) const
{
  int length;
  if(is_change_endian == true){
    Complex b[4];
    for(int i = 0; i < 2 * 4; i++)
      ((double*)b)[i] = ((double*)a)[i];
    yukawa::change_endian((double*)b, 2 * 4);
    length = fwrite(b, sizeof(Complex), 4, fp);
  }
  else{
    length = fwrite(a, sizeof(Complex), 4, fp);
  }
  if(length != 4){
    fprintf(stderr, "ERROR(PP::write): length = %d != 4\n", length);
    exit(1);
  }
}

/// change endian
inline PP& PP::change_endian()
{
  yukawa::change_endian((double*)a, 2 * 4);
  return *this;
}

//--------------------------------------------------------------------------
// helper functions
//--------------------------------------------------------------------------

inline PP operator+(PP lhs, const PP& rhs) { return lhs += rhs; }
inline PP operator-(PP lhs, const PP& rhs) { return lhs -= rhs; }
inline PP operator*(PP lhs, const PP& rhs) { return lhs *= rhs; }

inline PP operator-(PP pp)                 { return pp.unary_minus(); }

inline PP operator*(PP pp, double x)       { return pp *= x; }
inline PP operator*(double x, PP pp)       { return pp *= x; }
inline PP operator/(PP pp, double x)       { return pp /= x; }

inline PP operator*(PP pp, std::complex<double> z) { return pp *= z; }
inline PP operator*(std::complex<double> z, PP pp) { return pp *= z; }
inline PP operator/(std::complex<double> z, PP pp) { return pp /= z; }

inline PP rotate(PP pp, int ith) { return pp.rotate(ith); }
inline PP reflectioN(PP pp)      { return pp.reflection(); }

End_Name_Space_Yukawa

#endif
