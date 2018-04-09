//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup gauge
 * @brief   Header file for complex matrix class (arbitrary rank) & SU3 class
 * @author  Takaya Miyamoto
 * @since   Sat Dec 19 14:39:00 JST 2015
 */
//--------------------------------------------------------------------------

#ifndef COMPLEX_MATRIX_H
#define COMPLEX_MATRIX_H

#include <complex>
#include <cstddef>

#ifndef COMP_IMAG
#define COMP_IMAG std::complex<double>(0.0,1.0)
#endif
#ifndef COMP_ZERO
#define COMP_ZERO std::complex<double>(0.0,0.0)
#endif

//--------------------------------------------------------------------------
/**
 * @brief The class for complex matrix (arbitrary dimension)
 */
//--------------------------------------------------------------------------
class cmatrix
{
private:
   std::complex<double> *mat;
   size_t m_Ndim;
    
protected:
   
public:
//================== For inner index ==================//
   size_t ij(const size_t i, const size_t j)
   {
      return j + m_Ndim * i;
   }
   const size_t ij(const size_t i, const size_t j) const
   {
      return j + m_Ndim * i;
   }
//================== For writing ==================//
   std::complex<double>& operator()(const size_t i, const size_t j)
   {
      return mat[ij(i,j)];
   }
   std::complex<double>& operator()(const size_t i)
   {
      return mat[i];
   }
//================== For reading ==================//
   const std::complex<double>& operator()(const size_t i, const size_t j) const
   {
      return mat[ij(i,j)];
   }
   const std::complex<double>& operator()(const size_t i) const
   {
      return mat[i];
   }
//================== Constructor ==================//
   cmatrix()
   {
      mat = NULL;
   }
   cmatrix(const size_t Ndim)
   {
      mat = NULL;
      init(Ndim);
   }
   cmatrix(const cmatrix &obj)
   {
      mat = NULL;
      init(obj.m_Ndim);
      
      for (size_t i=0; i<m_Ndim*m_Ndim; i++) (*this)(i) = obj(i);
   }
   ~cmatrix()   // Destructor
   {
      if (mat != NULL) delete [] mat;
   }
//================== For initialize ==================//
   void init(const size_t Ndim)
   {
      if (mat != NULL) delete [] mat;
      m_Ndim = Ndim;
      mat = new std::complex<double>[m_Ndim*m_Ndim];
   }
   void init()
   {
      if (mat != NULL)
      {
         delete [] mat;
         mat = new std::complex<double>[m_Ndim*m_Ndim];
      }
      else
      {
         printf(" @ WARNING in cmatrix.init() : Must be set matrix dimension.\n");
      }
   }
//================== Operator define ==================//
   cmatrix& operator=(const cmatrix &rhs)
   {
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         (*this)(i) = rhs(i);
      
      return *this;
   }
   cmatrix& operator=(const std::complex<double> &rhs)
   {
      for (size_t i=0; i<m_Ndim; i++) for (size_t j=0; j<m_Ndim; j++)
         i==j ? (*this)(i,j) = rhs : (*this)(i,j) = COMP_ZERO;
      
      return *this;
   }
   cmatrix& operator+=(const cmatrix &rhs)
   {
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         (*this)(i) += rhs(i);
      
      return *this;
   }
   cmatrix& operator+=(const std::complex<double> &rhs)
   {
      for (size_t i=0; i<m_Ndim; i++)
         (*this)(i,i) += rhs;
      
      return *this;
   }
   cmatrix& operator-=(const cmatrix &rhs)
   {
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         (*this)(i) -= rhs(i);
      
      return *this;
   }
   cmatrix& operator-=(const std::complex<double> &rhs)
   {
      for (size_t i=0; i<m_Ndim; i++)
         (*this)(i,i) -= rhs;
      
      return *this;
   }
   cmatrix& operator*=(const cmatrix &rhs)
   {
      cmatrix tmp(m_Ndim);
      for (size_t i=0; i<m_Ndim; i++) for (size_t j=0; j<m_Ndim; j++)
         for (size_t k=0; k<m_Ndim; k++)
            tmp(i,j) += (*this)(i,k) * rhs(k,j);
      
      return *this = tmp;
   }
   cmatrix& operator*=(const std::complex<double> &rhs)
   {
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         (*this)(i) *= rhs;
      
      return *this;
   }
   cmatrix& operator/=(const std::complex<double> &rhs)
   {
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         (*this)(i) /= rhs;
      
      return *this;
   }
//================== Operator helper ==================//
   friend cmatrix operator+(cmatrix lhs, const cmatrix &rhs)
   {
      return lhs += rhs;
   }
   friend cmatrix operator+(cmatrix lhs, const std::complex<double> &rhs)
   {
      return lhs += rhs;
   }
   friend cmatrix operator+(const std::complex<double> &lhs, cmatrix rhs)
   {
      return rhs += lhs;
   }
   friend cmatrix operator-(cmatrix lhs, const cmatrix &rhs)
   {
      return lhs -= rhs;
   }
   friend cmatrix operator-(cmatrix lhs, const std::complex<double> &rhs)
   {
      return lhs -= rhs;
   }
   friend cmatrix operator-(const std::complex<double> &lhs, cmatrix rhs)
   {
      rhs *= std::complex<double>(-1.0);
      return rhs += lhs;
   }
   friend cmatrix operator*(cmatrix lhs, const cmatrix &rhs)
   {
      return lhs *= rhs;
   }
   friend cmatrix operator*(cmatrix lhs, const std::complex<double> &rhs)
   {
      return lhs *= rhs;
   }
   friend cmatrix operator*(const std::complex<double> &lhs, cmatrix rhs)
   {
      return rhs *= lhs;
   }
   friend cmatrix operator/(cmatrix lhs, const std::complex<double> &rhs)
   {
      return lhs /= rhs;
   }
//================== Several function ==================//
   size_t dimension()
   {
      return m_Ndim;
   }
   cmatrix unit()
   {
      cmatrix ret(m_Ndim);
      for (size_t i=0; i<m_Ndim; i++)
         ret(i,i) = std::complex<double>(1.0,0.0);
      
      return ret;
   }
   cmatrix dagger()
   {
      cmatrix ret(m_Ndim);
      for (size_t i=0; i<m_Ndim; i++) for (size_t j=0; j<m_Ndim; j++)
         ret(j,i) = conj((*this)(i,j));
      
      return ret;
   }
   std::complex<double> trace()
   {
      std::complex<double> ret;
      for (size_t i=0; i<m_Ndim; i++)
         ret += (*this)(i,i);
      
      return ret;
   }
   std::complex<double> sum_all_elements()
   {
      std::complex<double> ret;
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         ret += (*this)(i);
      
      return ret;
   }
   double sum_all_abs_elements()
   {
      double ret = 0.0;
      for (size_t i=0; i<m_Ndim*m_Ndim; i++)
         ret += abs((*this)(i));
      
      return ret;
   }
   cmatrix pow_mat(const int Npower)
   {
      cmatrix ret(*this);
      for (int n=1; n<Npower; n++)
         ret *= (*this);
      
      return ret;
   }
   cmatrix exp_mat(const int Dim_truncate)   // mat_ret = exp( mat_org )
   {
      cmatrix ret = (   std::complex<double>(1.0)
                     + (*this)/std::complex<double>(Dim_truncate) );
      
      for (int loop=Dim_truncate-1; loop>0; loop--)
         ret = (   std::complex<double>(1.0)
                + (*this)*ret/std::complex<double>(loop) );
      
      return ret;
   }
   void print()
   {
      for (size_t i=0; i<m_Ndim; i++)
      {
         for (size_t j=0; j<m_Ndim; j++)
            printf(" (%lf,%lf)", (*this)(i,j).real(), (*this)(i,j).imag());
         printf("\n");
      }
      printf("\n");
   }
};

//--------------------------------------------------------------------------
/**
 * @brief The class for complex matrix (For SU3)
 */
//--------------------------------------------------------------------------
class SU3matrix
{
private:
   std::complex<double> mat[9];
   
protected:
   
public:
//================== For inner index ==================//
   int ij(const int i, const int j)
   {
      return j + 3 * i;
   }
   const int ij(const int i, const int j) const
   {
      return j + 3 * i;
   }
//================== For writing ==================//
   std::complex<double>& operator()(const int i, const int j)
   {
      return mat[ij(i,j)];
   }
   std::complex<double>& operator()(const int i)
   {
      return mat[i];
   }
//================== For reading ==================//
   const std::complex<double>& operator()(const int i, const int j) const
   {
      return mat[ij(i,j)];
   }
   const std::complex<double>& operator()(const int i) const
   {
      return mat[i];
   }
//================== Constructor ==================//
   SU3matrix() {}
   SU3matrix(const SU3matrix &obj)
   {
      for (int i=0; i<9; i++) (*this)(i) = obj(i);
   }
   SU3matrix(const std::complex<double> &valu)
   {
      for (int i=0; i<3; i++) (*this)(i,i) = valu;
   }
   ~SU3matrix() {}   // Destructor
//================== For initialize ==================//
   void init()
   {
      for (int i=0; i<9; i++) (*this)(i) = COMP_ZERO;
   }
//================== Operator define ==================//
   SU3matrix& operator=(const SU3matrix &rhs)
   {
      for (int i=0; i<9; i++)
         (*this)(i) = rhs(i);
      
      return *this;
   }
   SU3matrix& operator=(const std::complex<double> &rhs)
   {
      for (int i=0; i<3; i++) for (int j=0; j<3; j++)
         i==j ? (*this)(i,j) = rhs : (*this)(i,j) = COMP_ZERO;
      
      return *this;
   }
   SU3matrix& operator+=(const SU3matrix &rhs)
   {
      for (int i=0; i<9; i++)
         (*this)(i) += rhs(i);
      
      return *this;
   }
   SU3matrix& operator+=(const std::complex<double> &rhs)
   {
      for (int i=0; i<3; i++)
         (*this)(i,i) += rhs;
      
      return *this;
   }
   SU3matrix& operator-=(const SU3matrix &rhs)
   {
      for (int i=0; i<9; i++)
         (*this)(i) -= rhs(i);
      
      return *this;
   }
   SU3matrix& operator-=(const std::complex<double> &rhs)
   {
      for (int i=0; i<3; i++)
         (*this)(i,i) -= rhs;
      
      return *this;
   }
   SU3matrix& operator*=(const SU3matrix &rhs)
   {
      SU3matrix tmp;
      for (int i=0; i<3; i++) for (int j=0; j<3; j++)
         for (int k=0; k<3; k++)
            tmp(i,j) += (*this)(i,k) * rhs(k,j);
      
      return *this = tmp;
   }
   SU3matrix& operator*=(const std::complex<double> &rhs)
   {
      for (int i=0; i<9; i++)
         (*this)(i) *= rhs;
      
      return *this;
   }
   SU3matrix& operator/=(const std::complex<double> &rhs)
   {
      for (int i=0; i<9; i++)
         (*this)(i) /= rhs;
      
      return *this;
   }
//================== Operator helper ==================//
   friend SU3matrix operator+(SU3matrix lhs, const SU3matrix &rhs)
   {
      return lhs += rhs;
   }
   friend SU3matrix operator+(SU3matrix lhs, const std::complex<double> &rhs)
   {
      return lhs += rhs;
   }
   friend SU3matrix operator+(const std::complex<double> &lhs, SU3matrix rhs)
   {
      return rhs += lhs;
   }
   friend SU3matrix operator-(SU3matrix lhs, const SU3matrix &rhs)
   {
      return lhs -= rhs;
   }
   friend SU3matrix operator-(SU3matrix lhs, const std::complex<double> &rhs)
   {
      return lhs -= rhs;
   }
   friend SU3matrix operator-(const std::complex<double> &lhs, SU3matrix rhs)
   {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend SU3matrix operator*(SU3matrix lhs, const SU3matrix &rhs)
   {
      return lhs *= rhs;
   }
   friend SU3matrix operator*(SU3matrix lhs, const std::complex<double> &rhs)
   {
      return lhs *= rhs;
   }
   friend SU3matrix operator*(const std::complex<double> &lhs, SU3matrix rhs)
   {
      return rhs *= lhs;
   }
   friend SU3matrix operator/(SU3matrix lhs, const std::complex<double> &rhs)
   {
      return lhs /= rhs;
   }
//================== Several function ==================//
   int dimension()
   {
      return 3;
   }
   static SU3matrix unit()
   {
      SU3matrix ret;
      for (int i=0; i<3; i++)
         ret(i,i) = std::complex<double>(1.0,0.0);
      
      return ret;
   }
   SU3matrix dagger()
   {
      SU3matrix ret;
      for (int i=0; i<3; i++) for (int j=0; j<3; j++)
         ret(j,i) = conj((*this)(i,j));
      
      return ret;
   }
   std::complex<double> trace()
   {
      return (*this)(0,0) + (*this)(1,1) + (*this)(2,2);
   }
   std::complex<double> det()
   {
      return  (  (*this)(0,0)*(*this)(1,1)*(*this)(2,2)
               + (*this)(1,0)*(*this)(2,1)*(*this)(0,2)
               + (*this)(2,0)*(*this)(0,1)*(*this)(1,2)
               - (*this)(0,2)*(*this)(1,1)*(*this)(2,0)
               - (*this)(0,1)*(*this)(1,0)*(*this)(2,2)
               - (*this)(0,0)*(*this)(2,1)*(*this)(1,2) );
   }
   SU3matrix inverce()
   {
      SU3matrix tmp;
      tmp(0,0) = (*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
      tmp(0,1) = (*this)(0,2)*(*this)(2,1)-(*this)(0,1)*(*this)(2,2);
      tmp(0,2) = (*this)(0,1)*(*this)(1,2)-(*this)(0,2)*(*this)(1,1);
      tmp(1,0) = (*this)(1,2)*(*this)(2,0)-(*this)(1,0)*(*this)(2,2);
      tmp(1,1) = (*this)(0,0)*(*this)(2,2)-(*this)(0,2)*(*this)(2,0);
      tmp(1,2) = (*this)(0,2)*(*this)(1,0)-(*this)(0,0)*(*this)(1,2);
      tmp(2,0) = (*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0);
      tmp(2,1) = (*this)(0,1)*(*this)(2,0)-(*this)(0,0)*(*this)(2,1);
      tmp(2,2) = (*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0);
      
      return tmp/(*this).det();
   }
   std::complex<double> sum_all_elements()
   {
      std::complex<double> ret;
      for (int i=0; i<9; i++)
         ret += (*this)(i);
      
      return ret;
   }
   double sum_all_abs_elements()
   {
      double ret = 0.0;
      for (int i=0; i<9; i++)
         ret += abs((*this)(i));
      
      return ret;
   }
   SU3matrix pow_mat(const int Npower)
   {
      SU3matrix ret(*this);
      for (int n=1; n<Npower; n++)
         ret *= (*this);
      
      return ret;
   }
   SU3matrix exp_mat(const int Dim_truncate)   // mat_ret = exp( mat_org )
   {
      SU3matrix ret = 1.0 + (*this)/std::complex<double>(Dim_truncate);
      
      for (int loop=Dim_truncate-1; loop>0; loop--)
         ret = 1.0 + (*this)*ret/std::complex<double>(loop);
      
      return ret;
   }
   SU3matrix exp_su3()   // mat_ret = exp( i * mat_org ),   (Exact value)
   {
      double c0     =   (*this).det().real();
      double c1     = (((*this)*(*this)).trace().real()) / 2.0;
      double c0_max = 2.0 * pow(c1/3.0, 1.5);
      double theta  = acos(c0/c0_max);
      double www    = sqrt(c1)     * sin(theta/3.0);
      double uuu    = sqrt(c1/3.0) * cos(theta/3.0);
      double q1     =  2.0 * uuu;
      double q2     = -uuu + www;
      double q3     = -uuu - www;
      double qdet   = q2*q3*q3+q3*q1*q1+q1*q2*q2-q2*q1*q1-q3*q2*q2-q1*q3*q3;
      
      std::complex<double> f0 = (  (q2*q3*q3-q3*q2*q2)*exp(COMP_IMAG*q1)
                                 + (q1*q1*q3-q1*q3*q3)*exp(COMP_IMAG*q2)
                                 + (q1*q2*q2-q2*q1*q1)*exp(COMP_IMAG*q3) ) / qdet;
      std::complex<double> f1 = (  (q2*q2   -q3*q3   )*exp(COMP_IMAG*q1)
                                 + (q3*q3   -q1*q1   )*exp(COMP_IMAG*q2)
                                 + (q1*q1   -q2*q2   )*exp(COMP_IMAG*q3) ) / qdet;
      std::complex<double> f2 = (  (q3      -q2      )*exp(COMP_IMAG*q1)
                                 + (q1      -q3      )*exp(COMP_IMAG*q2)
                                 + (q2      -q1      )*exp(COMP_IMAG*q3) ) / qdet;
      
      return f0 + (*this)*f1 + ((*this)*(*this))*f2;
   }
   void print()
   {
      for (int i=0; i<3; i++)
      {
         for (int j=0; j<3; j++)
            printf(" (%lf,%lf)", (*this)(i,j).real(), (*this)(i,j).imag());
         printf("\n");
      }
      printf("\n");
   }
};

#endif
