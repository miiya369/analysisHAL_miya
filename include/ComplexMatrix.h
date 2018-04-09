//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Matrix
 * @brief   Header file of class for arbitrary rank complex matrix & SU3 class
 * @author  Takaya Miyamoto
 * @since   Tue Nov  8 05:10:39 JST 2016
 */
//--------------------------------------------------------------------------

#ifndef COMPLEX_MATRIX_H
#define COMPLEX_MATRIX_H

#include "MatrixTemplate.h"
#include <complex>

typedef std::complex<double> cdouble;

#ifndef COMP_IMAG
#define COMP_IMAG cdouble(0.0,1.0)
#endif
#ifndef COMP_ZERO
#define COMP_ZERO cdouble(0.0,0.0)
#endif

//--------------------------------------------------------------------------
/**
 * @brief The class for the arbitrary dimension complex matrix
 */
//--------------------------------------------------------------------------
class cmatrix : public MATRIX_TEMPLATE_BASE<cdouble>
{
public:
   typedef MATRIX_TEMPLATE_BASE<cdouble> base;
   typedef cmatrix              myclass;
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//================== Constructor ==================//
   cmatrix() {}
   cmatrix(const myclass& other) : base(other) {}
   cmatrix(const    base& other) : base(other) {}
   cmatrix(const int Ndim)       : base(Ndim)  {}
   ~cmatrix() {}
//================== Operator definition ==================//
   myclass& operator =(const double &rhs) {
      for (int i=0; i<m_Ndim; i++) for (int j=0; j<m_Ndim; j++)
         i==j ? (*this)(i,j) = rhs : (*this)(i,j) = 0.0;
      return *this;
   }
   myclass& operator+=(const double &rhs) {
      for (int i=0; i<m_Ndim;        i++) (*this)(i,i) += rhs;
      return *this;
   }
   myclass& operator-=(const double &rhs) {
      for (int i=0; i<m_Ndim;        i++) (*this)(i,i) -= rhs;
      return *this;
   }
   myclass& operator*=(const double &rhs) {
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i)   *= rhs;
      return *this;
   }
   myclass& operator/=(const double &rhs) {
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i)   /= rhs;
      return *this;
   }
//================== Operator helper ==================//
   friend myclass operator+(myclass lhs, const myclass &rhs) {
      return lhs += rhs;
   }
   friend myclass operator-(myclass lhs, const myclass &rhs) {
      return lhs -= rhs;
   }
   friend myclass operator*(myclass lhs, const myclass &rhs) {
      return lhs *= rhs;
   }
   friend myclass operator+(myclass lhs, const cdouble &rhs) {
      return lhs += rhs;
   }
   friend myclass operator-(myclass lhs, const cdouble &rhs) {
      return lhs -= rhs;
   }
   friend myclass operator*(myclass lhs, const cdouble &rhs) {
      return lhs *= rhs;
   }
   friend myclass operator/(myclass lhs, const cdouble &rhs) {
      return lhs /= rhs;
   }
   friend myclass operator+(myclass lhs, const  double &rhs) {
      return lhs += rhs;
   }
   friend myclass operator-(myclass lhs, const  double &rhs) {
      return lhs -= rhs;
   }
   friend myclass operator*(myclass lhs, const  double &rhs) {
      return lhs *= rhs;
   }
   friend myclass operator/(myclass lhs, const  double &rhs) {
      return lhs /= rhs;
   }
   friend myclass operator+(const cdouble &lhs, myclass rhs) {
      return rhs += lhs;
   }
   friend myclass operator-(const cdouble &lhs, myclass rhs) {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend myclass operator*(const cdouble &lhs, myclass rhs) {
      return rhs *= lhs;
   }
   friend myclass operator+(const  double &lhs, myclass rhs) {
      return rhs += lhs;
   }
   friend myclass operator-(const  double &lhs, myclass rhs) {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend myclass operator*(const  double &lhs, myclass rhs) {
      return rhs *= lhs;
   }
//================== Several function ==================//
   myclass T() {
      return base::T();
   }
   myclass get_unit() {
      return base::get_unit();
   }
   myclass inverce() {
      return base::inverce();
   }
   myclass pow_mat(const int Npower) {
      return base::pow_mat(Npower);
   }
   myclass dagger() {
      myclass ret(m_Ndim);
      for (int i=0; i<m_Ndim; i++) for (int j=0; j<m_Ndim; j++)
         ret(j,i) = conj((*this)(i,j));
      return ret;
   }
   double sum_all_abs_elements() {
      double ret = 0.0;
      for (int i=0; i<m_Ndim*m_Ndim; i++) ret += abs((*this)(i));
      return ret;
   }
   myclass exp_mat(const int Dim_truncate) {   // mat_ret = exp( mat_org )
      myclass ret = (1.0 + (*this)/double(Dim_truncate));
      for (int loop=Dim_truncate-1; loop>0; loop--)
         ret = (1.0 + ret*(*this)/double(loop));
      return ret;
   }
   void print() {
      for (   int i=0; i<m_Ndim; i++) {
         for (int j=0; j<m_Ndim; j++)
            printf("(%12lf,%12lf) ", (*this)(i,j).real(), (*this)(i,j).imag());
         printf("\n");
      }
      printf("\n");
   }
   void print_e() {
      for (   int i=0; i<m_Ndim; i++) {
         for (int j=0; j<m_Ndim; j++)
            printf("(%15e,%15e) ", (*this)(i,j).real(), (*this)(i,j).imag());
         printf("\n");
      }
      printf("\n");
   }
};

//--------------------------------------------------------------------------
/**
 * @brief The generators for su3 matrix
 */
//--------------------------------------------------------------------------
namespace generator {
   cmatrix su2(const int);
   cmatrix su3(const int);
}

//--------------------------------------------------------------------------
/**
 * @brief The class for su3 matrix
 */
//--------------------------------------------------------------------------
class SU3matrix : public cmatrix
{
public:
   typedef cmatrix   base;
   typedef SU3matrix myclass;
//================== Constructor ==================//
   SU3matrix()                     : base(3)     {}
   SU3matrix(const myclass& other) : base(other) {}
   SU3matrix(const    base& other) : base(other) {
      if (other.Ndim() != 3) ERROR_COMMENTS("Cast to SU3mat from cmatrix "
                                            "with Ndim != 3 is not allowed");
   }
   SU3matrix(const  cdouble& valu) : base(3) {
      for (int i=0; i<3; i++) for (int j=0; j<3; j++)
         i==j ? (*this)(i,j) = valu : (*this)(i,j) = 0.0;
   }
   SU3matrix(const   double& valu) : base(3) {
      for (int i=0; i<3; i++) for (int j=0; j<3; j++)
         i==j ? (*this)(i,j) = valu : (*this)(i,j) = 0.0;
   }
   ~SU3matrix() {}
//================== Several function ==================//
   base exp_su3() {   // mat_ret = exp( i * mat_org ),   (Exact value)
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
      
      cdouble f0 = ((q2*q3*q3-q3*q2*q2)*exp(COMP_IMAG*q1) +
                    (q1*q1*q3-q1*q3*q3)*exp(COMP_IMAG*q2) +
                    (q1*q2*q2-q2*q1*q1)*exp(COMP_IMAG*q3)) / qdet;
      cdouble f1 = ((q2*q2   -q3*q3   )*exp(COMP_IMAG*q1) +
                    (q3*q3   -q1*q1   )*exp(COMP_IMAG*q2) +
                    (q1*q1   -q2*q2   )*exp(COMP_IMAG*q3)) / qdet;
      cdouble f2 = ((q3      -q2      )*exp(COMP_IMAG*q1) +
                    (q1      -q3      )*exp(COMP_IMAG*q2) +
                    (q2      -q1      )*exp(COMP_IMAG*q3)) / qdet;
      
      return f0 + (*this)*f1 + ((*this)*(*this))*f2;
   }
   //! return Exp(i * Coeffs[i] * Gell-Mann matrix[i+1]), i=0-7
   base make_SU3(const double Coeffs[8]) {
      myclass ret(3); ret.set_zeros();
      for (int i=0; i<8; i++) ret += (generator::su3(i+1) * Coeffs[i]);
      return ret.exp_su3();
   }
};

#endif
