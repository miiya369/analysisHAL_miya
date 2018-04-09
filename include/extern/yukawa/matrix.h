//--------------------------------------------------------------------------
/**
 * @file
 * @brief   matrix of variable dimension
 * @ingroup General_Util
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Fri May 21 22:48:40 JST 2010
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_MATRIX_H
#define IS_INCLUDED_MATRIX_H

#include <stdio.h>
#include <cstdlib>
#include <complex>

#include <yukawa/yukawa.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

//--------------------------------------------------------------------------
/**
 * @brief matrix of variable dimension
 * @since Wed Feb 24 02:49:42 JST 2010
 */
//--------------------------------------------------------------------------

template<class T>
class matrix {
  typedef T         Scalar;
  typedef matrix<T> Matrix;

  int    _size;     //< dimension of the matrix
  Scalar *a;
public:
  ~matrix() { if(a) delete[] a; }
  matrix();
  matrix(int n);
  matrix(const Scalar ary[],int dim);
  matrix(const Matrix& other);

  Matrix& operator=(const Matrix& other);

  int     size()                  const { return _size; }
  const
  Scalar& operator()(int i,int j) const { return a[range(i,j)]; }
  Scalar& operator()(int i,int j)       { return a[range(i,j)]; }

  Matrix& operator+=(const Matrix& other);
  Matrix& operator-=(const Matrix& other);
  Matrix& operator*=(const Matrix& other);

  Matrix& operator*=(Scalar z);
  Matrix& operator/=(Scalar z);

  Scalar  tr()        const;
  Matrix  transpose() const;
  Matrix  dagger()    const;

private:
  int range(int i,int j) const;
  void check_size(const Matrix& other) const;
};

template<class T> inline
int matrix<T>::range(int i,int j) const
{
  if(!NDEBUG){
    if(i < 0 || size() <= i || j < 0 || size() <= j){
      fprintf(stderr, "ERROR(range) in matrix<T>\n");
      exit(1);
    }
  }
#if 1
  return j + size()*i;
#else
  return i + size()*j;
#endif
}

template<class T> inline
void matrix<T>::check_size(const Matrix& other) const
{
  if(!NDEBUG){
    if(size() != other.size()){
      fprintf(stderr, "ERROR(size) in matrix<T>\n");
      exit(1);
    }
  }
}

/// the default constructor
template <class T> inline
matrix<T>::matrix()
  : _size(0), a(0)
{}

/// constructor 1
template <class T> inline
matrix<T>::matrix(int n)
  : _size(n), a(0)
{
  if(_size != 0){
    a = new Scalar[_size*_size];
    for(int i = 0; i < _size*_size; i++) a[i] = 0.0;
  }
}

/// constructor 2
template <class T> inline
matrix<T>::matrix(const T ary[],int dim)
  : _size(dim), a(0)
{
  if(_size != 0){
    a = new Scalar[_size*_size];
    for(int i = 0; i < _size*_size; i++) a[i] = ary[i];
  }
}

/// the copy constructor
template <class T> inline
matrix<T>::matrix(const matrix<T>& other)
  : _size(other._size), a(0)
{
  if(_size != 0){
    a = new Scalar[_size*_size];
    for(int i = 0; i < _size*_size; i++) a[i] = other.a[i];
  }
}

/// the assignment operator
template <class T> inline
matrix<T>& matrix<T>::operator=(const matrix<T>& other)
{
  _size = other._size;
  if(a) delete[] a;
  a = 0;

  if(_size != 0){
    a = new Scalar[_size*_size];
    for(int i = 0; i < _size*_size; i++) a[i] = other.a[i];
  }
  return *this;
}

/// += operation
template <class T> inline
matrix<T>& matrix<T>::operator+=(const matrix<T>& other)
{
  check_size(other);
  for(int i = 0; i < size()*size(); i++) a[i] += other.a[i];
  return *this;
}

/// -= operation
template <class T> inline
matrix<T>& matrix<T>::operator-=(const matrix<T>& other)
{
  check_size(other);
  for(int i = 0; i < size()*size(); i++) a[i] -= other.a[i];
  return *this;
}

/// *= operation (matrix multiplication)
template <class T> inline
matrix<T>& matrix<T>::operator*=(const matrix<T>& other)
{
  check_size(other);
  Matrix self(*this);
  for(  int i = 0; i < size(); i++)
    for(int j = 0; j < size(); j++){
      Scalar sum = 0.0;
      for(int k = 0; k < size(); k++)
	sum += self(i,k) * other(k,j);

      (*this)(i,j) = sum;
    }
  return *this;
}

/// *= operation (scalar multiplication)
template <class T> inline
matrix<T>& matrix<T>::operator*=(T z)
{
  for(int i = 0; i < size()*size(); i++) a[i] *= z;
  return *this;
}

/// /= operation (scalar division)
template <class T> inline
matrix<T>& matrix<T>::operator/=(T z)
{
  for(int i = 0; i < size()*size(); i++) a[i] /= z;
  return *this;
}

/// trace
template <class T> inline
T matrix<T>::tr() const
{
  T sum = T(0.0);
  for(int i = 0; i < size(); i++)
    sum += (*this)(i,i);

  return sum;
}

/// transpose
template <class T> inline
matrix<T> matrix<T>::transpose() const
{
  Matrix ret(size());
  for(  int i = 0; i < size(); i++)
    for(int j = 0; j < size(); j++)
      ret(i,j) = (*this)(j,i);

  return ret;
}

/// dagger
template <class T> inline
matrix<T> matrix<T>::dagger() const
{
  Matrix ret(size());
  for(  int i = 0; i < size(); i++)
    for(int j = 0; j < size(); j++)
      ret(i,j) = conj((*this)(j,i));

  return ret;
}


//--------------------------------------------------------------------------
// Helper functions
//--------------------------------------------------------------------------

/// matrix addition
template <class T> inline
matrix<T> operator+(matrix<T> lhs,const matrix<T>& rhs)
{ return lhs += rhs; }

/// matrix subtraction
template <class T> inline
matrix<T> operator-(matrix<T> lhs,const matrix<T>& rhs)
{ return lhs -= rhs; }

/// matrix multiplication
template <class T> inline
matrix<T> operator*(matrix<T> lhs,const matrix<T>& rhs)
{ return lhs *= rhs; }

/// scalar multiplication 1
template <class T> inline
matrix<T> operator*(matrix<T> mat,T z)
{ return mat *= z; }

/// scalar multiplication 2
template <class T> inline
matrix<T> operator*(T z,matrix<T> mat)
{ return mat *= z; }

/// scalar division
template <class T> inline
matrix<T> operator/(matrix<T> mat,T z)
{ return mat /= z; }

/// trace
template <class T> inline
T tr(const matrix<T>& mat) { return mat.tr(); }

/// transpose
template <class T> inline
matrix<T> transpose(const matrix<T>& mat) { return mat.transpose(); }

/// dagger
template <class T> inline
matrix<T> dagger(const matrix<T>& mat) { return mat.dagger(); }

/// operator <<
template <class T> inline
std::ostream& operator<<(std::ostream& os,const matrix<T>& mat)
{
  for(  int i=0; i<mat.size(); i++){
    for(int j=0; j<mat.size(); j++){
      os << mat(i,j) << "\t";
    }
    os << "\n";
  }
  return os;
}

//---------------------------------------------------------------------
/**
 * @brief returns the Kronecker product of two matrices
 *
 * Let @f$a@f$  and @f$b@f$ be n by  n and m by  m matrices. Kronecker
 * product  of  @f$a@f$  and  @f$b@f$  is  a  (n  x  m)  by  (n  x  m)
 * matrix,denoted by @f$a\otimes b@f$ defined by
 *
 * @f[
 *   (a \otimes b)_{i,j} = a_{i_1,j_1} b_{i_2,j_2}
 * @f]
 * with @f$i = i_1 + n i_2@f$ and @f$j  = j_1 + n j_2$.
 */
//---------------------------------------------------------------------

template <class T>
inline matrix<T>
Kronecker_product(const matrix<T>& mat1,const matrix<T>& mat2)
{
  matrix<T> ret(mat1.size() * mat2.size());

  for(  int i1 = 0; i1 < mat1.size(); i1++)
    for(int j1 = 0; j1 < mat1.size(); j1++)
      for(  int i2 = 0; i2 < mat2.size(); i2++)
	for(int j2 = 0; j2 < mat2.size(); j2++){
	  ret(i1 + mat1.size()*i2,j1 + mat1.size()*j2)
	    = mat1(i1,j1)*mat2(i2,j2);
	}

  return ret;
}

//--------------------------------------------------------------------------
// special operations for matrix<Complex>
//--------------------------------------------------------------------------

/// +=: matrix<Complex>, matrix<double>
inline matrix<std::complex<double> >&
operator+=(matrix<std::complex<double> >& lhs,
	   const matrix<double>& rhs)
{
  for(  int j = 0; j < lhs.size(); j++)
    for(int i = 0; i < lhs.size(); i++)
      lhs(i,j) += rhs(i,j);
  return lhs;
}

/// -=: matrix<Complex>, matrix<double>
inline matrix<std::complex<double> >&
operator-=(matrix<std::complex<double> >& lhs,
	   const matrix<double>& rhs)
{
  for(  int j = 0; j < lhs.size(); j++)
    for(int i = 0; i < rhs.size(); i++)
      lhs(i,j) -= rhs(i,j);
  return lhs;
}

/// *= matrix<Complex>, double
inline matrix<std::complex<double> >&
operator*=(matrix<std::complex<double> >& mat,
	   double x)
{
  for(  int i = 0; i < mat.size(); i++)
    for(int j = 0; j < mat.size(); j++)
      mat(i,j) *= x;
  return mat;
}

/// /= matrix<Complex>, double
inline matrix<std::complex<double> >&
operator/=(matrix<std::complex<double> >& mat,
	   double x)
{
  for(  int i = 0; i < mat.size(); i++)
    for(int j = 0; j < mat.size(); j++)
      mat(i,j) /= x;
  return mat;
}

/// addition: matrix<Complex> + matrix<double>
inline matrix<std::complex<double> >
operator+(matrix<std::complex<double> > lhs,const matrix<double>& rhs)
{
  for(  int j = 0; j < lhs.size(); j++)
    for(int i = 0; i < lhs.size(); i++)
      lhs(i,j) += rhs(i,j);
  return lhs;
}

/// addition: matrix<double> + matrix<Complex>
inline matrix<std::complex<double> >
operator+(const matrix<double>& lhs,matrix<std::complex<double> > rhs)
{
  for(  int j = 0; j < lhs.size(); j++)
    for(int i = 0; i < lhs.size(); i++)
      rhs(i,j) += rhs(i,j);
  return rhs;
}

/// subtraction: matrix<Complex> - matrix<double>
inline matrix<std::complex<double> >
operator-(matrix<std::complex<double> > lhs,const matrix<double>& rhs)
{
  for(  int j = 0; j < lhs.size(); j++)
    for(int i = 0; i < lhs.size(); i++)
      lhs(i,j) -= rhs(i,j);
  return lhs;
}

/// addition: matrix<double> - matrix<Complex>
inline matrix<std::complex<double> >
operator-(const matrix<double>& lhs,matrix<std::complex<double> > rhs)
{
  for(  int j = 0; j < lhs.size(); j++)
    for(int i = 0; i < lhs.size(); i++)
      rhs(i,j) = lhs(i,j) - rhs(i,j);
  return rhs;
}

//--------------------------------------------------------------------------
// another spetial operations for matrix<Complex>
//--------------------------------------------------------------------------

inline matrix<std::complex<double> >
operator*(double x,matrix<std::complex<double> > mat)
{ return mat *= x; }

inline matrix<std::complex<double> >
operator*(matrix<std::complex<double> > mat,double x)
{ return mat *= x; }

inline matrix<std::complex<double> >
operator/(matrix<std::complex<double> > mat,double x)
{ return mat /= x; }

// eigen values and eigen vectors

extern void jacobi(const matrix<double>& mat,
		   double d[],
		   double v[],
		   int& nrot);

extern void jacobi(const matrix<std::complex<double> >& mat,
		   std::complex<double> d[],
		   std::complex<double> v[],
		   int& nrot);

// conversion

extern
matrix<std::complex<double> > ftoc(const matrix<double>& mat);

// I/O

inline
void fprint(FILE* fp, const char* fmt, const matrix<double>& mat)
{
  for(  int i = 0; i < mat.size(); i++){
    for(int j = 0; j < mat.size(); j++){
      fprintf(fp, fmt, mat(i,j));
    }
    fprintf(fp, "\n");
  }
}

inline
void fprint(FILE* fp, const char* fmt, const matrix<std::complex<double> >& mat)
{
  for(  int i = 0; i < mat.size(); i++){
    for(int j = 0; j < mat.size(); j++){
      fprintf(fp, fmt, real(mat(i,j)), imag(mat(i,j)));
    }
    fprintf(fp, "\n");
  }  
}

End_Name_Space_Yukawa

#endif
