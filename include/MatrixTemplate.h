//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Matrix
 * @brief   Header file of template base class for arbitrary dimension matrix
 * @author  Takaya Miyamoto
 * @since   Wed Nov  1 16:44:07 JST 2017
 */
//--------------------------------------------------------------------------

#ifndef MATRIX_TEMPLATE_BASE_H
#define MATRIX_TEMPLATE_BASE_H

#include <stdio.h>
#include <stdlib.h>

#ifndef ERROR_COMMENTS
#define ERROR_COMMENTS(COMMENTS) { \
printf("\n\nERROR(%s, Line:%d): %s\n\n", __FILE__, __LINE__, COMMENTS); \
exit(1); \
}
#endif

//--------------------------------------------------------------------------
/**
 * @brief The template base class for arbitrary dimension matrix
 */
//--------------------------------------------------------------------------
template <class X> class MATRIX_TEMPLATE_BASE
{
protected:
   X*  m_mat;
   int m_Ndim;
   
//================== For inner index ==================//
   const int ij(const int i, const int j) const {
      return j + m_Ndim * i;
   }
public:
   typedef MATRIX_TEMPLATE_BASE myclass;
//================== For writing ==================//
   X& operator()(const int i, const int j) {
      if (i >= m_Ndim || j >= m_Ndim) ERROR_COMMENTS("Index Overflow");
      return m_mat[ij(i,j)];
   }
   X& operator()(const int index) {
      if (index >= m_Ndim*m_Ndim)     ERROR_COMMENTS("Index Overflow");
      return m_mat[index];
   }
//================== For reading ==================//
   const X& operator()(const int i, const int j) const {
      if (i >= m_Ndim || j >= m_Ndim) ERROR_COMMENTS("Index Overflow");
      return m_mat[ij(i,j)];
   }
   const X& operator()(const int index) const {
      if (index >= m_Ndim*m_Ndim)     ERROR_COMMENTS("Index Overflow");
      return m_mat[index];
   }
//================== Constructor ==================//
   MATRIX_TEMPLATE_BASE() {
      m_mat  = NULL;
      m_Ndim = 0;
   }
   MATRIX_TEMPLATE_BASE(const myclass& other) {
      m_mat   = NULL;
      m_Ndim  = 0;
      (*this) = other;
   }
   MATRIX_TEMPLATE_BASE(const int Ndim) {
      m_mat  = NULL;
      m_Ndim = 0;
      init(Ndim);
   }
   ~MATRIX_TEMPLATE_BASE() {
      if (m_mat != NULL) delete [] m_mat;
   }
//================== For initialize ==================//
   void init(const int Ndim) {
      if (m_Ndim != Ndim) {
         if (m_mat != NULL) {
            delete [] m_mat;
            m_mat = NULL;
         }
         m_Ndim = Ndim;
      }
      if (m_mat == NULL)
         m_mat = new X[m_Ndim * m_Ndim];
   }
//================== Operator definition ==================//
   myclass& operator =(const myclass &rhs) {
      init(rhs.Ndim());
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i)  = rhs(i);
      return *this;
   }
   myclass& operator+=(const myclass &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i) += rhs(i);
      return *this;
   }
   myclass& operator-=(const myclass &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i) -= rhs(i);
      return *this;
   }
   myclass& operator*=(const myclass &rhs) {
      if ((*this).Ndim() != rhs.Ndim()) ERROR_COMMENTS("Different dimensions.");
      myclass tmp(*this);
      for (int i=0; i<m_Ndim; i++) for (int j=0; j<m_Ndim; j++) {
         (*this)(i,j) = 0.0;
         for (int k=0; k<m_Ndim; k++) (*this)(i,j) += tmp(i,k) * rhs(k,j);
      }
      return *this;
   }
   myclass& operator =(const X &rhs) {
      for (int i=0; i<m_Ndim; i++) for (int j=0; j<m_Ndim; j++)
         i==j ? (*this)(i,j) = rhs : (*this)(i,j) = 0.0;
      return *this;
   }
   myclass& operator+=(const X &rhs) {
      for (int i=0; i<m_Ndim;        i++) (*this)(i,i) += rhs;
      return *this;
   }
   myclass& operator-=(const X &rhs) {
      for (int i=0; i<m_Ndim;        i++) (*this)(i,i) -= rhs;
      return *this;
   }
   myclass& operator*=(const X &rhs) {
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i)   *= rhs;
      return *this;
   }
   myclass& operator/=(const X &rhs) {
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
   friend myclass operator+(const X &lhs, myclass rhs) {
      return rhs += lhs;
   }
   friend myclass operator-(const X &lhs, myclass rhs) {
      int l_Ndim = rhs.Ndim();
      for (int i=0; i<l_Ndim*l_Ndim; i++) rhs(i) *= (-1.0);
      return rhs += lhs;
   }
   friend myclass operator*(const X &lhs, myclass rhs) {
      return rhs *= lhs;
   }
//================== Several functions ==================//
   const int Ndim() const {
      return m_Ndim;
   }
   void set_zeros() {
      for (int i=0; i<m_Ndim*m_Ndim; i++) (*this)(i) = 0.0;
   }
   void set_unit() {
      for (int i=0; i<m_Ndim; i++) for (int j=0; j<m_Ndim; j++)
         i==j ? (*this)(i,j) = 1.0 : (*this)(i,j) = 0.0;
   }
   myclass T() {
      myclass ret(m_Ndim);
      for (int i=0; i<m_Ndim; i++) for (int j=0; j<m_Ndim; j++)
         ret(j,i) = (*this)(i,j);
      return ret;
   }
   myclass get_unit() {
      myclass ret(*this);
      ret.set_unit();
      return ret;
   }
   X trace() {
      X ret((*this)(0,0));
      for (int i=1; i<m_Ndim; i++) ret += (*this)(i,i);
      return ret;
   }
   X det() {
      if      (m_Ndim == 1)
         return  (*this)(0,0);
      else if (m_Ndim == 2)
         return ((*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0));
      else if (m_Ndim == 3)
         return ((*this)(0,0)*(*this)(1,1)*(*this)(2,2) +
                 (*this)(1,0)*(*this)(2,1)*(*this)(0,2) +
                 (*this)(2,0)*(*this)(0,1)*(*this)(1,2) -
                 (*this)(0,2)*(*this)(1,1)*(*this)(2,0) -
                 (*this)(0,1)*(*this)(1,0)*(*this)(2,2) -
                 (*this)(0,0)*(*this)(2,1)*(*this)(1,2));
      else
         ERROR_COMMENTS("myclass.det() at "
                        "Ndim > 3 is not implemented yet.\n");
   }
   myclass inverce() {
      myclass tmp((*this).Ndim());
      if      (m_Ndim == 1) {
         tmp(0,0) =  (*this)(0,0); tmp(0,0) = 1.0;
      }
      else if (m_Ndim == 2) {
         tmp(0,0) = (*this)(1,1);
         tmp(0,1) = (*this)(0,1) * (-1.0);
         tmp(1,0) = (*this)(1,0) * (-1.0);
         tmp(1,1) = (*this)(0,0);
      }
      else if (m_Ndim == 3) {
         tmp(0,0) = (*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
         tmp(0,1) = (*this)(0,2)*(*this)(2,1)-(*this)(0,1)*(*this)(2,2);
         tmp(0,2) = (*this)(0,1)*(*this)(1,2)-(*this)(0,2)*(*this)(1,1);
         tmp(1,0) = (*this)(1,2)*(*this)(2,0)-(*this)(1,0)*(*this)(2,2);
         tmp(1,1) = (*this)(0,0)*(*this)(2,2)-(*this)(0,2)*(*this)(2,0);
         tmp(1,2) = (*this)(0,2)*(*this)(1,0)-(*this)(0,0)*(*this)(1,2);
         tmp(2,0) = (*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0);
         tmp(2,1) = (*this)(0,1)*(*this)(2,0)-(*this)(0,0)*(*this)(2,1);
         tmp(2,2) = (*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0);
      }
      else
         ERROR_COMMENTS("myclass.det() at "
                        "Ndim > 3 is not implemented yet.\n");
      tmp /= (*this).det();
      return tmp;
   }
   X sum_all_elements() {
      X ret((*this)(0));
      for (int i=1; i<m_Ndim*m_Ndim; i++) ret += (*this)(i);
      return ret;
   }
   myclass pow_mat(const int Npower) {
      myclass ret(*this);
      for (int n=1; n<Npower; n++) ret *= (*this);
      return ret;
   }
};

#endif
