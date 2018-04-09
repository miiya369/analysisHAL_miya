//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Matrix
 * @brief   Header file for class of complex field matrix
 * @author  Takaya Miyamoto
 * @since   Tue Nov  8 05:10:39 JST 2016
 */
//--------------------------------------------------------------------------

#ifndef COMPLEX_FIELD_MATRIX_H
#define COMPLEX_FIELD_MATRIX_H

#include "ComplexMatrix.h"

//--------------------------------------------------------------------------
/**
 * @brief The template class of complex field matrix
 */
//--------------------------------------------------------------------------
template <class X> class ComplexFieldMatrix : public MATRIX_TEMPLATE_BASE<X>
{
public:
   typedef MATRIX_TEMPLATE_BASE<X> base;
   typedef ComplexFieldMatrix      myclass;
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//================== Constructor ==================//
   ComplexFieldMatrix() {}
   ComplexFieldMatrix(const myclass& other) : base(other) {}
   ComplexFieldMatrix(const    base& other) : base(other) {}
   ComplexFieldMatrix(const int Ndim)       : base(Ndim)  {}
   ~ComplexFieldMatrix() {}
//================== Operator definition ==================//
   myclass& operator =(const cmatrix &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i)  = rhs(i);
      return *this;
   }
   myclass& operator+=(const cmatrix &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i) += rhs(i);
      return *this;
   }
   myclass& operator-=(const cmatrix &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i) -= rhs(i);
      return *this;
   }
   myclass& operator*=(const cmatrix &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      myclass tmp(*this);
      for (int i=0; i<base::m_Ndim; i++) for (int j=0; j<base::m_Ndim; j++) {
         (*this)(i,j) = 0.0;
         for (int k=0; k<base::m_Ndim; k++) (*this)(i,j) += tmp(i,k) * rhs(k,j);
      }
      return *this;
   }
   myclass& operator =(const cdouble &rhs) {
      for (int i=0; i<base::m_Ndim; i++) for (int j=0; j<base::m_Ndim; j++)
         i==j ? (*this)(i,j) = rhs : (*this)(i,j) = 0.0;
      return *this;
   }
   myclass& operator+=(const cdouble &rhs) {
      for (int i=0; i<base::m_Ndim;              i++) (*this)(i,i) += rhs;
      return *this;
   }
   myclass& operator-=(const cdouble &rhs) {
      for (int i=0; i<base::m_Ndim;              i++) (*this)(i,i) -= rhs;
      return *this;
   }
   myclass& operator*=(const cdouble &rhs) {
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i)   *= rhs;
      return *this;
   }
   myclass& operator/=(const cdouble &rhs) {
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i)   /= rhs;
      return *this;
   }
   myclass& operator =(const double &rhs) {
      for (int i=0; i<base::m_Ndim; i++) for (int j=0; j<base::m_Ndim; j++)
         i==j ? (*this)(i,j) = rhs : (*this)(i,j) = 0.0;
      return *this;
   }
   myclass& operator+=(const double &rhs) {
      for (int i=0; i<base::m_Ndim;              i++) (*this)(i,i) += rhs;
      return *this;
   }
   myclass& operator-=(const double &rhs) {
      for (int i=0; i<base::m_Ndim;              i++) (*this)(i,i) -= rhs;
      return *this;
   }
   myclass& operator*=(const double &rhs) {
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i)   *= rhs;
      return *this;
   }
   myclass& operator/=(const double &rhs) {
      for (int i=0; i<base::m_Ndim*base::m_Ndim; i++) (*this)(i)   /= rhs;
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
   friend myclass operator+(myclass lhs, const cmatrix &rhs) {
      return lhs += rhs;
   }
   friend myclass operator-(myclass lhs, const cmatrix &rhs) {
      return lhs -= rhs;
   }
   friend myclass operator*(myclass lhs, const cmatrix &rhs) {
      return lhs *= rhs;
   }
   friend myclass operator+(myclass lhs, const X &rhs) {
      return lhs += rhs;
   }
   friend myclass operator-(myclass lhs, const X &rhs) {
      return lhs -= rhs;
   }
   friend myclass operator*(myclass lhs, const X &rhs) {
      return lhs *= rhs;
   }
   friend myclass operator/(myclass lhs, const X &rhs) {
      return lhs /= rhs;
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
   friend myclass operator+(const cmatrix &lhs, myclass rhs) {
      return rhs += lhs;
   }
   friend myclass operator-(const cmatrix &lhs, myclass rhs) {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend myclass operator*(const cmatrix &lhs, myclass rhs) {
      if (lhs.Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      myclass ret(rhs);
      int l_Ndim = ret.Ndim();
      for (int i=0; i<l_Ndim; i++) for (int j=0; j<l_Ndim; j++) {
         ret(i,j) = 0.0;
         for (int k=0; k<l_Ndim; k++) ret(i,j) += lhs(i,k) * rhs(k,j);
      }
      return ret;
   }
   friend myclass operator+(const X &lhs, myclass rhs) {
      return rhs += lhs;
   }
   friend myclass operator-(const X &lhs, myclass rhs) {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend myclass operator*(const X &lhs, myclass rhs) {
      return rhs *= lhs;
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
      myclass ret(base::m_Ndim);
      for (int i=0; i<base::m_Ndim; i++) for (int j=0; j<base::m_Ndim; j++)
         ret(j,i) = (*this)(i,j).conj();
      return ret;
   }
};

#endif
