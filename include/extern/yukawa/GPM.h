//---------------------------------------------------------------------
/**
 * @file
 * @brief   A generalization of permutation matrices (GPM)
 * @ingroup General_Util
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Thu May 22 23:13:06 JST 2008
 *
 * Each row of a permutation matrix has only one entry, whose value is
 * one.  It  is easy to implement  it with a table.  For instance, for
 * matrix<br>
 *
 * @f[ 
 * A \equiv
 *
 * \left[
 * \begin{array}{cccc}
 * 0 & 1 & 0 & 0 \\
 * 1 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 1 \\
 * 0 & 0 & 1 & 0
 * \end{array}
 * \right]
 *@f]
 * 
 * @f[
 * \begin{array}{c|cccc}
 *    & 0 & 1 & 2 & 3 \\
 *  \hline
 *  I & 1 & 0 & 3 & 2
 * \end{array}
 * @f]
 *
 * @f$ A_{ij} = \delta_{j,I(i)}. @f$
 *
 * We generalize it so that each row has only one entry, but its value
 * is not necessarily one. In this case, a matrix is represented by an
 * extended table. For instance, for matrix<br>
 * @f[ 
 * A \equiv
 *
 * \left[
 * \begin{array}{cccc}
 * 0 & a & 0 & 0 \\
 * b & 0 & 0 & 0 \\
 * 0 & 0 & 0 & c \\
 * 0 & 0 & d & 0
 * \end{array}
 * \right]
 *@f]
 * 
 * @f[
 * \begin{array}{c|cccc}
 *    & 0 & 1 & 2 & 3 \\
 *  \hline
 *  I & 1 & 0 & 3 & 2 \\
 *  Z & a & b & c & d 
 * \end{array}
 * @f]
 *
 * @f$ A_{ij} = Z(i) \delta_{j,I(i)}. @f$
 */
//---------------------------------------------------------------------

#ifndef IS_INCLUDED_GPM_H
#define IS_INCLUDED_GPM_H

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <unistd.h>

#include <yukawa/yukawa.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

inline double dagger(double x) { return x; }
inline double conj(  double x) { return x; }

inline int    dagger(int x) { return x; }
inline int    conj(  int x) { return x; }

//---------------------------------------------------------------------
/**
 * @brief A generalization of permutation matrices (GPM)
 *
 *        @f[ A_{ij} = Z(i) \delta_{j,I(i)} @f]
 */
//---------------------------------------------------------------------

template<class T>
struct GPM {
  typedef T Scalar;

  int  size;
  int* I;
  T*   Z;

  ~GPM() { if(I) delete[] I; if(Z) delete[] Z; }
  GPM() : size(0), I(0), Z(0) {}
  GPM(int i0,int i1,               T t0,T t1);
  GPM(int i0,int i1,int i2,        T t0,T t1,T t2);
  GPM(int i0,int i1,int i2,int i3, T t0,T t1,T t2,T t3);
  GPM(int dim);
  GPM(int dim, int i[], const T t[]);
  GPM(const GPM<T>& other)
    : size(other.size), I(new int[size]), Z(new T[size])
  { for(int i = 0; i < size; i++){ I[i] = other.I[i]; Z[i] = other.Z[i]; } }

  GPM<T>& operator=(const GPM<T>& other){
    size = other.size;
    if(I) delete[] I; if(Z) delete[] Z;
    I = new int[size]; Z = new T[size];
    for(int i = 0; i < size; i++){ I[i] = other.I[i]; Z[i] = other.Z[i]; }
    return *this;
  }

  void resize(int n){
    if(I) delete[] I; if(Z) delete[] Z;
    size = n;
    I = new int[size];
    Z = new T[size];
  }

  T operator()(int i,int j) const;

  GPM<T>& operator*=(const GPM<T>& other){
    check_size(other);
    for(int i = 0; i < size; i++){
      Z[i] = Z[i] * other.Z[I[i]];
      I[i] = other.I[I[i]];
    }
    return *this;
  }

  void transform(const T src[],T dst[]) const;

  void check_size(const GPM<T>& other) const{
    if(!NDEBUG)
      if(size != other.size){
	fprintf(stderr,
		"ERROR(GPM): inconsistent size\n");
	exit(1);
      }
  }
};

//---------------------------------------------------------------------
/**
 * @brief a constructor for an object with dimension two
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T>::GPM(int i0,int i1, T t0,T t1)
  : size(2), I(new int[2]), Z(new T[2])
{
  I[0]=i0; I[1]=i1;
  Z[0]=t0; Z[1]=t1;
}

//---------------------------------------------------------------------
/**
 * @brief a constructor for an object with dimension three
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T>::GPM(int i0,int i1,int i2, T t0,T t1,T t2)
  : size(3), I(new int[3]), Z(new T[3])
{
  I[0]=i0; I[1]=i1; I[2]=i2;
  Z[0]=t0; Z[1]=t1; Z[2]=t2;
}

//---------------------------------------------------------------------
/**
 * @brief a constructor for an object with dimension four
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T>::GPM(int i0, int i1, int i2, int i3,
	    T   t0, T   t1, T   t2, T   t3)
  : size(4), I(new int[4]), Z(new T[4])
{
  I[0]=i0; I[1]=i1; I[2]=i2; I[3]=i3;
  Z[0]=t0; Z[1]=t1; Z[2]=t2; Z[3]=t3;
}

//---------------------------------------------------------------------
/**
 * @brief a constructor for an object with general dimension
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T>::GPM(int dim, int i[], const T t[])
  : size(dim), I(new int[dim]), Z(new T[dim])
{
  for(int j = 0; j < size; j++){
    I[j] = i[j];
    Z[j] = t[j];
  }
}

//---------------------------------------------------------------------
/**
 * @brief
 */
//---------------------------------------------------------------------

template<class T> inline
T GPM<T>::operator()(int i,int j) const
{
  if(j != I[i]) return T();
  return Z[i];
}

//---------------------------------------------------------------------
/**
 * @brief action on vectors of class T
 */
//---------------------------------------------------------------------

template<class T> inline
void GPM<T>::transform(const T src[],T dst[]) const
{
  for(int i = 0; i < size; i++) dst[i] = Z[i]*src[I[i]];
}

//---------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------

//---------------------------------------------------------------------
/**
 * @brief print the content in matrix form
 */
//---------------------------------------------------------------------

template<class T> inline
std::ostream& operator<<(std::ostream& os,const GPM<T>& obj)
{
  for(int i = 0; i < obj.size; i++){
    for(int j = 0; j < obj.size; j++){
      if(j == obj.I[i])
	os << obj.Z[i] << "\t";
      else
	os << 0 << "\t";
    }
    os << "\n";
  }
  return os;
}

//---------------------------------------------------------------------
/**
 * @brief matrix multiplication
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T>
operator*(GPM<T> lhs, const GPM<T>& rhs)
{ return lhs *= rhs; }

//---------------------------------------------------------------------
/**
 * @brief Kroneker product of matrices
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T>
Kronecker_product(const GPM<T>& lhs, const GPM<T>& rhs)
{
  GPM<T> ret;
  ret.resize(lhs.size * rhs.size);

  for(  int i = 0; i < lhs.size; i++)
    for(int j = 0; j < rhs.size; j++){
      ret.Z[i + lhs.size*j] = lhs.Z[i] * rhs.Z[j];
      ret.I[i + lhs.size*j] = lhs.I[i] + lhs.size*rhs.I[j];
    }
  return ret;
}

//---------------------------------------------------------------------
/**
 * @brief transposed matrix
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T> transpose(const GPM<T>& gpm)
{
  GPM<T> ret(gpm.size);
  for(int i = 0; i < ret.size; i++){
    ret.Z[gpm.I[i]] = gpm.Z[i];
    ret.I[gpm.I[i]] = i;
  }
  return ret;
}

//---------------------------------------------------------------------
/**
 * @brief dagger
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T> dagger(const GPM<T>& gpm)
{
  GPM<T> ret;
  ret.resize(gpm.size);

  for(int i = 0; i < ret.size; i++){
    ret.Z[gpm.I[i]] = conj(gpm.Z[i]);
    ret.I[gpm.I[i]] = i;
  }
  return ret;
}

//---------------------------------------------------------------------
/**
 * @brief inverse
 */
//---------------------------------------------------------------------

template<class T> inline
GPM<T> inverse(const GPM<T>& gpm)
{
  GPM<T> ret;
  ret.resize(gpm.size);

  for(int i = 0; i < ret.size; i++){
    ret.Z[gpm.I[i]] = T(1)/gpm.Z[i];
    ret.I[gpm.I[i]] = i;
  }
  return ret;
}

End_Name_Space_Yukawa

#endif


