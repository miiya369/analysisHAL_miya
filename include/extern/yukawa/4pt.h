//--------------------------------------------------------------------------
/**
 * @file
 * @brief
 * @ingroup Four_Point_Function
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Fri May 21 18:52:36 JST 2010
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup Four_Point_Function Four Point Function
 * @ingroup  General_Util
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_4PT_H
#define IS_INCLUDED_4PT_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/GPM.h>
#include <yukawa/cubic_group.h>

#include <yukawa.CORE/PP.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

struct   four_point_base {
  typedef std::complex<double> Complex;
  enum BC {PBC=0, APBC=1};

  ~four_point_base(){}
  static void initialize() { cubic_group::initialize(); }
};

//--------------------------------------------------------------------------
/**
 * @brief four point correlator of hadron fields
 *
 *
 */
//--------------------------------------------------------------------------

template <class T>
struct   four_point : public four_point_base {

  int Xsites, Ysites, Zsites, XYZsites;
  BC  Xbc,    Ybc,    Zbc;
  
  T* a;

  ~four_point();
  four_point();
  four_point(int xsites,  int ysites,  int zsites,
	     BC  xbc=PBC, BC  ybc=PBC, BC zbc=PBC);
  four_point(const four_point<T>& other);
  
  four_point<T>& operator=(const four_point<T>& other);

  four_point<T>& operator+=(const four_point<T>& other);
  four_point<T>& operator-=(const four_point<T>& other);

  four_point<T>& unary_minus();

  four_point<T>& operator*=(double x);
  four_point<T>& operator*=(Complex z);
  four_point<T>& operator*=(const four_point<T>& other);

  four_point<T>& operator/=(double x);
  four_point<T>& operator/=(Complex z);

  const
  T& operator()(int ix,int iy,int iz, int& isign) const;
  T& operator()(int ix,int iy,int iz, int& isign);

  const
  T& operator[](int ixyz) const { return a[range(ixyz)]; }
  T& operator[](int ixyz)       { return a[range(ixyz)]; }

  four_point<T>& rotate(int ith, bool is_restricted_to_sink_xyz = false);
  four_point<T>& reflection(     bool is_restricted_to_sink_xyz = false);
  four_point<T>& charge_conjugation();

  void read( const char* fname, bool is_change_endian=false);
  void read( FILE* fp,          bool is_change_endian=false);

  void write(const char* fname, bool is_change_endian=false);
  void write(FILE* fp,          bool is_change_endian=false);

  four_point<PP> to_PP(Complex eta[]) const;

  // The following is related to N.Ishii, H.Nemura, and K.Murano
  four_point<T>& cps2standard();
private:
  int range(int ixyz) const;
  int range(int ix, int iy, int iz) const;

  void is_same_size(const four_point<T>& other) const;
};

/// range
template <class T> inline
int four_point<T>::range(int ixyz) const
{
  if(!NDEBUG){
    if(ixyz < 0 || XYZsites <= ixyz){
      fprintf(stderr,
	      "ERROR(four_point::range): ixyz/XYZsites=%d/%d\n",
	      ixyz,
	      XYZsites);
      //      exit(1);
      throw *this;
    }
  }
  return ixyz;
}

/// range 2
template <class T> inline
int four_point<T>::range(int ix, int iy, int iz) const
{
  if(!NDEBUG){
    if(ix < 0 || Xsites <= ix){
      fprintf(stderr,"ERROR(four_point::range): ix/Xsites=%d/%d\n",ix,Xsites);
      exit(1);
    }
    if(iy < 0 || Ysites <= iy){
      fprintf(stderr,"ERROR(four_point::range): iy/Ysites=%d/%d\n",iy,Ysites);
      exit(1);
    }
    if(iz < 0 || Zsites <= iz){
      fprintf(stderr,"ERROR(four_point::range): iz/Zsites=%d/%d\n",iz,Zsites);
      exit(1);
    }
  }
  return ix + Xsites*(iy + Ysites*iz);
}

/// size checker
template <class T> inline
void four_point<T>::is_same_size(const four_point<T>& other) const
{
  if(!NDEBUG){
    if(Xsites != other.Xsites){
      fprintf(stderr,
	      "ERROR(four_point::is_same_size): Xsites, %d != %d\n",
	      Xsites, other.Xsites);
      exit(1);
    }
    if(Ysites != other.Ysites){
      fprintf(stderr,
	      "ERROR(four_point::is_same_size): Ysites, %d != %d\n",
	      Ysites, other.Ysites);
      exit(1);
    }
    if(Zsites != other.Zsites){
      fprintf(stderr,
	      "ERROR(four_point::is_same_size): Zsites, %d != %d\n",
	      Zsites, other.Zsites);
      exit(1);
    }

    if(Xbc != other.Xbc){
      fprintf(stderr,
	      "ERROR(four_point::is_same_size): Xbc\n");
      exit(1);
    }
    if(Ybc != other.Ybc){
      fprintf(stderr,
	      "ERROR(four_point::is_same_size): Ybc\n");
      exit(1);
    }
    if(Zbc != other.Zbc){
      fprintf(stderr,
	      "ERROR(four_point::is_same_size): Zbc\n");
      exit(1);
    }
  }
}

/// destructor
template <class T> inline
four_point<T>::~four_point()
{ if(a) delete[] a; }

/// the default constructor
template <class T> inline
four_point<T>::four_point()
  : Xsites(0), Ysites(0), Zsites(0), XYZsites(0),
    Xbc(PBC),  Ybc(PBC),  Zbc(PBC),
    a(0)
{
  T::initialize();
}

/// a constructor
template <class T> inline
four_point<T>::four_point(int xsites, int ysites, int zsites,
			  BC  xbc,    BC  ybc,    BC  zbc)
  : Xsites(xsites), Ysites(ysites), Zsites(zsites),
    XYZsites(xsites * ysites * zsites),
    Xbc(xbc), Ybc(ybc), Zbc(zbc),
    a(new T[xsites * ysites * zsites])
{
  T::initialize();
}

/// the copy constructor
template <class T> inline
four_point<T>::four_point(const four_point<T>& other)
  : Xsites(other.Xsites), Ysites(other.Ysites), Zsites(other.Zsites),
    XYZsites(other.XYZsites),
    Xbc(other.Xbc), Ybc(other.Ybc), Zbc(other.Zbc),
    a(new T[Xsites * Ysites * Zsites])
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] = other.a[ixyz];
}

/// the assignment operator
template <class T> inline
four_point<T>& four_point<T>::operator=(const four_point<T>& other)
{
  if(a) delete[] a;
  Xsites = other.Xsites;
  Ysites = other.Ysites;
  Zsites = other.Zsites;

  XYZsites = other.XYZsites;

  Xbc = other.Xbc;
  Ybc = other.Ybc;
  Zbc = other.Zbc;

  a = new T[XYZsites];

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] = other.a[ixyz];

  return *this;
}

/// += operation
template <class T> inline
four_point<T>& four_point<T>::operator+=(const four_point<T>& other)
{
  is_same_size(other);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] += other.a[ixyz];
  return *this;
}

/// -= operation
template <class T> inline
four_point<T>& four_point<T>::operator-=(const four_point<T>& other)
{
  is_same_size(other);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] -= other.a[ixyz];
  return *this;
}

/// unary_minus
template <class T> inline
four_point<T>& four_point<T>::unary_minus()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] = -a[ixyz];
  return *this;
}

/// scalar multiplication (double)
template <class T> inline
four_point<T>& four_point<T>::operator*=(double x)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] *= x;
  return *this;
}

/// scalar multiplication (Complex)
template <class T> inline
four_point<T>& four_point<T>::operator*=(Complex x)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] *= x;
  return *this;
}

/// component-wise multiplication (double)
template <class T> inline
four_point<T>& four_point<T>::operator*=(const four_point<T>& other)
{
  is_same_size(other);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] *= other.a[ixyz];
  return *this;
}

/// scalar division (double)
template <class T> inline
four_point<T>& four_point<T>::operator/=(double x)
{
  double rec = 1.0/x;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] *= rec;
  return *this;
}

/// scalar division (Complex)
template <class T> inline
four_point<T>& four_point<T>::operator/=(Complex x)
{
  Complex rec = 1.0/x;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz] *= rec;
  return *this;
}

/// reference
template <class T> inline
const
T& four_point<T>::operator()(int ix, int iy, int iz, int& isign) const
{
  ix += 1024*Xsites;
  iy += 1024*Ysites;
  iz += 1024*Zsites;

  isign = 1;
  if(Xbc == APBC){
    if((ix / Xsites) % 2 == 1) isign = -isign;
  }
  if(Ybc == APBC){
    if((iy / Ysites) % 2 == 1) isign = -isign;
  }
  if(Zbc == APBC){
    if((iz / Zsites) % 2 == 1) isign = -isign;
  }

  ix = ix % Xsites;
  iy = iy % Ysites;
  iz = iz % Zsites;

  int ixyz = ix + Xsites*(iy + Ysites*iz);

  return a[ixyz];
}

/// reference 2
template <class T> inline
T& four_point<T>::operator()(int ix, int iy, int iz, int& isign)
{
  ix += 1024*Xsites;
  iy += 1024*Ysites;
  iz += 1024*Zsites;

  isign = 1;
  if(Xbc == APBC){
    if((ix / Xsites) % 2 == 1) isign = -isign;
  }
  if(Ybc == APBC){
    if((iy / Ysites) % 2 == 1) isign = -isign;
  }
  if(Zbc == APBC){
    if((iz / Zsites) % 2 == 1) isign = -isign;
  }

  ix = ix % Xsites;
  iy = iy % Ysites;
  iz = iz % Zsites;

  int ixyz = ix + Xsites*(iy + Ysites*iz);

  return a[ixyz];
}


/// rotation
template <class T> inline
four_point<T>& four_point<T>::rotate(int ith, bool is_restricted_to_xyz)
{
  using namespace cubic_group;

  four_point<T> self(*this);

  GPM<int> ginv(inverse(cube[ith]));

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
    int x[4], ginv_x[4];
    if(is_restricted_to_xyz == false){
      // full rotation
#ifdef _OPENMP
#pragma omp for
      for(int ixyz = 0; ixyz < XYZsites; ixyz++){
	int tmp = ixyz;
	x[0] = tmp % Xsites; tmp /= Xsites;
	x[1] = tmp % Ysites; tmp /= Ysites;
	x[2] = tmp;
#else
	for(    x[2] = 0; x[2] < Zsites; x[2]++)
	  for(  x[1] = 0; x[1] < Ysites; x[1]++)
	    for(x[0] = 0; x[0] < Xsites; x[0]++){
	      int ixyz = x[0] + Xsites*(x[1] + Ysites*x[2]);
#endif
#if 0
	    }
#endif
	ginv.transform(x,ginv_x);
	int sign=1;
	T* dst = &a[ixyz];
	T* src = &self(ginv_x[0],ginv_x[1],ginv_x[2], sign);
	
	if(sign==1){
	  *dst = (*src).rotate(ith);
	}
	else{
	  *dst = -(*src).rotate(ith);
	}
      }
    }
    else{
      // partial rotation with respect to xyz (sink)
#ifdef _OPENMP
#pragma omp for
      for(int ixyz = 0; ixyz < XYZsites; ixyz++){
	int tmp = ixyz;
	x[0] = tmp % Xsites; tmp /= Xsites;
	x[1] = tmp % Ysites; tmp /= Ysites;
	x[2] = tmp;
#else
	for(    x[2] = 0; x[2] < Zsites; x[2]++)
	  for(  x[1] = 0; x[1] < Ysites; x[1]++)
	    for(x[0] = 0; x[0] < Xsites; x[0]++){
	      int ixyz = x[0] + Xsites*(x[1] + Ysites*x[2]);
#endif
#if 0
	    }
#endif
	ginv.transform(x,ginv_x);
	int sign=1;
	T* dst = &a[ixyz];
	T* src = &self(ginv_x[0],ginv_x[1],ginv_x[2], sign);
	
	if(sign==1){
	  *dst = (*src);
	}
	else{
	  *dst = -(*src);
	}
      }
    }
#ifdef _OPENMP
  }
#endif
  return *this;
}

/// spatial reflection
template <class T> inline
four_point<T>& four_point<T>::reflection(bool is_restricted_to_xyz)
{
  four_point<T> self(*this);

  if(is_restricted_to_xyz == false){
    // full reflection
#ifdef _OPENMP
#pragma omp parallel for
    for(int ixyz = 0; ixyz < XYZsites; ixyz++){
      int tmp = ixyz;
      int ix = tmp % Xsites; tmp /= Xsites;
      int iy = tmp % Ysites; tmp /= Ysites;
      int iz = tmp;
#else
    for(    int iz = 0; iz < Zsites; iz++)
      for(  int iy = 0; iy < Ysites; iy++)
	for(int ix = 0; ix < Xsites; ix++){
#endif
	  int isign;
	  T*  src = &self(-ix,-iy,-iz,isign);
	  
	  if(isign == 1){
	    a[range(ix,iy,iz)] = (*src).reflection();
	  }
	  else{
	    a[range(ix,iy,iz)] = -(*src).reflection();
	  }
	}
  }
  else{
    // partial reflection with respect to xyz (sink)
#ifdef _OPENMP
#pragma omp parallel for
    for(int ixyz = 0; ixyz < XYZsites; ixyz++){
      int tmp = ixyz;
      int ix = tmp % Xsites; tmp /= Xsites;
      int iy = tmp % Ysites; tmp /= Ysites;
      int iz = tmp;
#else
    for(    int iz = 0; iz < Zsites; iz++)
      for(  int iy = 0; iy < Ysites; iy++)
	for(int ix = 0; ix < Xsites; ix++){
#endif
	  int isign;
	  T*  src = &self(-ix,-iy,-iz,isign);
	  
	  if(isign == 1){
	    a[range(ix,iy,iz)] = (*src);
	  }
	  else{
	    a[range(ix,iy,iz)] = -(*src);
	  }
	}
  }
  return *this;
}

/// charge_conjugation
template <class T> inline
four_point<T>& four_point<T>::charge_conjugation()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz].charge_conjugation();
  return *this;
}

/// to_PP (conversion to BS wave function)
template <class T> inline
four_point<PP> four_point<T>::to_PP(std::complex<double> eta[]) const
{
  four_point<PP> ret(Xsites, Ysites, Zsites, Xbc, Ybc, Zbc);
  for(int ixyz = 0; ixyz < XYZsites; ixyz++)
    ret[ixyz] = (*this)[ixyz].to_PP(eta);
  return ret;
}

/// cps2standard
template <class T> inline
four_point<T>& four_point<T>::cps2standard()
{
  for(int ixyz = 0; ixyz < XYZsites; ixyz++) a[ixyz].cps2standard();
  return *this;
}

//--------------------------------------------------------------------------
// helper functions
//--------------------------------------------------------------------------

/// addition
template <class T> inline
four_point<T> operator+(four_point<T> lhs,const four_point<T>& rhs)
{ return lhs += rhs; }

/// subtraction
template <class T> inline
four_point<T> operator-(four_point<T> lhs,const four_point<T>& rhs)
{ return lhs -= rhs; }

/// unary minus
template <class T> inline
four_point<T> operator-(four_point<T> four)
{ return four.unary_minus(); }



/// scalar multiplication (double)
template <class T> inline
four_point<T> operator*(four_point<T> four, double x)
{ return four *= x; }

/// scalar multiplication (double) 2
template <class T> inline
four_point<T> operator*(double x, four_point<T> four)
{ return four *= x; }

/// scalar multiplication (Complex)
template <class T> inline
four_point<T> operator*(four_point<T> four, std::complex<double> x)
{ return four *= x; }

/// scalar multiplication (Complex) 2
template <class T> inline
four_point<T> operator*(std::complex<double> x, four_point<T> four)
{ return four *= x; }


/// scalar division (double)
template <class T> inline
four_point<T> operator/(four_point<T> four, double x)
{ return four /= x; }

/// scalar division (Complex)
template <class T> inline
four_point<T> operator/(four_point<T> four, std::complex<double> x)
{ return four /= x; }

/// rotation
template <class T> inline
four_point<T> rotate(four_point<T> four, int ith, bool is_restricted_to_xyz = false)
{ return four.rotate(ith, is_restricted_to_xyz); }

/// reflection
template <class T> inline
four_point<T> reflection(four_point<T> four, bool is_restricted_to_xyz = false)
{ return four.reflection(is_restricted_to_xyz); }

/// chage conjugation
template <class T> inline
four_point<T> chage_conjugation(four_point<T> four)
{ return four.charge_conjugation(); }

//--------------------------------------------------------------------------
// noise reductions
//--------------------------------------------------------------------------

/// noise reduction by spatial reflection
template <class T> inline
four_point<T> noise_reduction_reflection(four_point<T> four)
{
  four_point<T> tmp(four);
  tmp += four.reflection();
  //fprintf(stderr, "\tnoise reduction(reflection)\n");

  return tmp/2.0;
}

//--------------------------------------------------------------------------
/**
 * @brief noise reduction by rotation
 *
 * Efficiency is achieved by using the formula
 * @f[
 * \begin{array}{lcl}
 * \sum_{i=0}^{23} c_i
 * &=&
 * \left(1 + c_x + c_x^2 + c_x^3 + c_y + c_y^3\right)
 * \left(1 + c_z + c_z^2 + c_z^3\right)
 * \\
 * &=&
 * \left\{ (1 + c_x) (1 + c_x^2) + c_y + c_y^2 \right\}
 * (1+c_z) (1+c_z^2)
 * \end{array}
 * \]
 */
//--------------------------------------------------------------------------

template <class T> inline
four_point<T> noise_reduction_rotation(four_point<T> four)
#if 1
{
  // (1)
  // further improved efficiency
  //
  // ((1 + c_x) (1 + c_x^2) + c_y + c_y^3) (1 + c_z) (1 + c_z^2)
  //

  four_point<T> tmp(four);
  //fprintf(stderr, "\tnoise reduction(rotation):");

  // (1 + cz^2)
  tmp += rotate(four, 17); // cz^2
  four = tmp;

  // (1 + cz)
  tmp += rotate(four, 3);  // cz
  four = tmp;

  four_point<T> tmp2(four);

  // (1 + cx^2)
  tmp += rotate(four, 15); // cx^2
  four = tmp;

  // (1 + cx)
  tmp += rotate(four, 1);  // cx

  // +cy
  tmp += rotate(tmp2, 2);  // cy
  
  // +cy^3
  tmp += rotate(tmp2, 5);  // cy^3

  //fprintf(stderr, " DONE\n");

  return tmp/24.0;
}
#else
#if 1
{
  // (2)
  // improved efficiency by coset decomposition
  //
  // (1 + c_x + c_x^2 + c_x^3 + c_y + c_y^3) (1 + c_z + c_z^2 + c_z^3)
  //

  four_point<T> tmp(four);
  //fprintf(stderr, "\tnoise reduction(rotation):");

  tmp += rotate(four, 3);  // cz
  tmp += rotate(four, 6);  // cz^{-1}
  tmp += rotate(four, 17); // cz^2

  four = tmp;

  tmp += rotate(four, 1);  // cx
  tmp += rotate(four, 2);  // cy
  tmp += rotate(four, 4);  // cx^{-1}
  tmp += rotate(four, 5);  // cy^{-1}
  tmp += rotate(four, 15); // cx^2

  //fprintf(stderr, " DONE\n");

  return tmp/24.0;
}
#else
{
  // (3)
  // brute force version
  //
  // \sum_{i=0}^{23} c_i
  //

  four_point<T> tmp(four);
  fprintf(stderr, "\tnoise reduction(rotation):");
  for(int ith = 1; ith < 24; ith++){
    tmp += rotate(four, ith);
    fprintf(stderr, " %d", ith);
  }
  fprintf(stderr, "\n");

  return tmp/24.0;
}
#endif
#endif

//--------------------------------------------------------------------------
/**
 * @brief I/O(read)
 */
//--------------------------------------------------------------------------

template <class T> inline
void four_point<T>::read(const char* fname, bool is_change_endian)
{
  FILE* fp;
  bool  is_popen = false;

  if(strcmp(&fname[strlen(fname)-3],".gz")==0){
    is_popen = true;
  }
  else{
    is_popen = false;
  }

  if(is_popen == true){
    char str[2048];
    sprintf(str, "zcat %s", fname);
    if((fp = popen(str, "r"))==NULL){
      fprintf(stderr, "cannot popen '%s'\n", str);
      exit(1);
    }
  }
  else{
    if((fp = fopen(fname, "r"))==NULL){
      fprintf(stderr, "cannot fopen '%s'\n", fname);
      exit(1);
    }
  }
  read(fp, is_change_endian);
  if(is_popen == true){
    pclose(fp);
  }
  else{
    fclose(fp);
  }
}

//--------------------------------------------------------------------------
/**
 * @brief I/O(read 2)
 */
//--------------------------------------------------------------------------

template <class T>
inline void four_point<T>::read(FILE* fp, bool is_change_endian)
{
  int32_t l[4];
  fread(l, sizeof(int32_t), 4, fp);
  Xbc = Ybc = Zbc = PBC;
  Xsites = l[0];
  Ysites = l[1];
  Zsites = l[2];
  if(Xsites < 0){ Xbc = APBC; Xsites = -Xsites; }
  if(Ysites < 0){ Ybc = APBC; Ysites = -Ysites; }
  if(Zsites < 0){ Zbc = APBC; Zsites = -Zsites; }

  XYZsites = Xsites * Ysites * Zsites;

  if(a) delete[] a;
  a = new T[Xsites * Ysites * Zsites];

  for(int ixyz = 0; ixyz < XYZsites; ixyz++)
    (*this)[ixyz].read(fp,is_change_endian);
}

//--------------------------------------------------------------------------
/**
 * @brief I/O(write)
 */
//--------------------------------------------------------------------------

template <class T>
inline void four_point<T>::write(const char* fname, bool is_change_endian)
{
  FILE* fp;
  bool  is_popen = false;

  if(strcmp(&fname[strlen(fname)-3],".gz")==0){
    is_popen = true;
  }
  else{
    is_popen = false;
  }
  
  if(is_popen == true){
    char str[2048];
    sprintf(str, "gzip -c > %s", fname);
    if((fp = popen(str, "w"))==NULL){
      fprintf(stderr, "cannot popen '%s'\n", str);
      exit(1);
    }
  }
  else{
    if((fp = fopen(fname, "w"))==NULL){
      fprintf(stderr, "cannot fopen '%s'\n", fname);
      exit(1);
    }
  }
  write(fp, is_change_endian);
  if(is_popen == true){
    pclose(fp);
  }
  else{
    fclose(fp);
  }
}

//--------------------------------------------------------------------------
/**
 * @brief I/O (write 2)
 */
//--------------------------------------------------------------------------

template <class T>
inline void four_point<T>::write(FILE* fp, bool is_change_endian)
{
  int32_t l[4];
  if(Xbc==PBC) l[0] = Xsites; else l[0] = -Xsites;
  if(Ybc==PBC) l[1] = Ysites; else l[1] = -Ysites;
  if(Zbc==PBC) l[2] = Zsites; else l[2] = -Zsites;
  l[3] = 0;

  fwrite(l, sizeof(int32_t), 4, fp);
  for(int ixyz = 0; ixyz < XYZsites; ixyz++)
    (*this)[ixyz].write(fp, is_change_endian);
}

End_Name_Space_Yukawa

#endif

