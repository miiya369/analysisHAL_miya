//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Definitions for sub class of complex-field
 * @author  Takaya Miyamoto
 * @since   Mon Oct 30 18:50:56 JST 2017
 */
//--------------------------------------------------------------------------

#ifndef COMPLEX_FIELD_SUB_H
#define COMPLEX_FIELD_SUB_H

#include "ComplexField_Base.h"
#include "StatisticsTemplate.h"

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in time (only t index)
 * @brief complex double Field(t)
 * @brief t: time coodination
 */
//--------------------------------------------------------------------------
class ComplexField_T : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_T    myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//======================== Constructor & Destructor ======================//
   ComplexField_T() {}
   ComplexField_T(const myclass& other) : base(other) {}
   ComplexField_T(const base&    other) : base(1, 1, other.get_tSIZE(), 1) {
      (*this) = other;
   }
   ComplexField_T(const base& other, const int ia, const int ix, const int iy,
                  const int iz, const int ib)
   : base(1, 1, other.get_tSIZE(), 1) {
      for (int t=0; t<(*this).data_size(); t++)
         (*this)(t) = other(ia,ix,iy,iz,t,ib);
   }
   ComplexField_T(const base& other, const int ia, const int ixyz, const int ib)
   : base(1, 1, other.get_tSIZE(), 1) {
      for (int t=0; t<(*this).data_size(); t++)
         (*this)(t) = other(ia,ixyz,t,ib);
   }
   ComplexField_T(const int a_tSIZE) : base(1, 1, a_tSIZE, 1) {}
   
   ~ComplexField_T() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_tSIZE) { base::mem_alloc(1, 1, a_tSIZE, 1); }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_aSIZE() != 1 || rhs.get_xSIZE() != 1 || rhs.get_ySIZE() != 1 ||
          rhs.get_zSIZE() != 1 || rhs.get_bSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index a,x,y,z,b != 1 "
                        "cannot cast to ComplexField_T");
      mem_alloc(rhs.get_tSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_T; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_T; }
};

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in time
 * @brief complex double Field(a, t, b)
 * @brief a: inner degree of freedom
 * @brief t: time  coodination
 * @brief b: outer degree of freedom
 */
//--------------------------------------------------------------------------
class ComplexField_ATB : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_ATB  myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int a, const int t, const int b) {
      return m_field[aVtb(a,0,t,b)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int a, const int t, const int b) const {
      return m_field[aVtb(a,0,t,b)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_ATB() {}
   ComplexField_ATB(const myclass& other) : base(other) {}
   ComplexField_ATB(const base&    other) : base(other.get_aSIZE(), 1,
                                                 other.get_tSIZE(),
                                                 other.get_bSIZE()) {
      (*this) = other;
   }
   ComplexField_ATB(const base& other, const int ix, const int iy, const int iz)
   : base(other.get_aSIZE(), 1, other.get_tSIZE(), other.get_bSIZE()) {
      for (      int b=0; b<(*this).get_bSIZE(); b++)
         for (   int t=0; t<(*this).get_tSIZE(); t++)
            for (int a=0; a<(*this).get_aSIZE(); a++)
               (*this)(a,t,b) = other(a,ix,iy,iz,t,b);
   }
   ComplexField_ATB(const base& other, const int ixyz)
   : base(other.get_aSIZE(), 1, other.get_tSIZE(), other.get_bSIZE()) {
      for (      int b=0; b<(*this).get_bSIZE(); b++)
         for (   int t=0; t<(*this).get_tSIZE(); t++)
            for (int a=0; a<(*this).get_aSIZE(); a++)
               (*this)(a,t,b) = other(a,ixyz,t,b);
   }
   ComplexField_ATB(const int a_aSIZE, const int a_tSIZE, const int a_bSIZE)
   : base(a_aSIZE, 1, a_tSIZE, a_bSIZE) {}
   
   ~ComplexField_ATB() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_aSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, 1, a_tSIZE, a_bSIZE);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_xSIZE() != 1 || rhs.get_ySIZE() != 1 || rhs.get_zSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index x,y,z != 1 "
                        "cannot cast to ComplexField_ATB");
      mem_alloc(rhs.get_aSIZE(), rhs.get_tSIZE(), rhs.get_bSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_ATB; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_ATB; }
};

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in time (only at index)
 * @brief complex double Field(a, t)
 * @brief a: inner degree of freedom
 * @brief t: time  coodination
 */
//--------------------------------------------------------------------------
class ComplexField_AT : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_AT   myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int a, const int t) {
      return m_field[aVtb(a,0,t,0)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int a, const int t) const {
      return m_field[aVtb(a,0,t,0)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_AT() {}
   ComplexField_AT(const myclass& other) : base(other) {}
   ComplexField_AT(const base&    other) : base(other.get_aSIZE(), 1,
                                                other.get_tSIZE(), 1) {
      (*this) = other;
   }
   ComplexField_AT(const base& other, const int ix, const int iy, const int iz,
                   const int ib)
   : base(other.get_aSIZE(), 1, other.get_tSIZE(), 1) {
      for (   int t=0; t<(*this).get_tSIZE(); t++)
         for (int a=0; a<(*this).get_aSIZE(); a++)
            (*this)(a,t) = other(a,ix,iy,iz,t,ib);
   }
   ComplexField_AT(const base& other, const int ixyz, const int ib)
   : base(other.get_aSIZE(), 1, other.get_tSIZE(), 1) {
      for (   int t=0; t<(*this).get_tSIZE(); t++)
         for (int a=0; a<(*this).get_aSIZE(); a++)
            (*this)(a,t) = other(a,ixyz,t,ib);
   }
   ComplexField_AT(const int a_aSIZE, const int a_tSIZE)
   : base(a_aSIZE, 1, a_tSIZE, 1) {}
   
   ~ComplexField_AT() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_aSIZE, const int a_tSIZE) {
      base::mem_alloc(a_aSIZE, 1, a_tSIZE, 1);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_xSIZE() != 1 || rhs.get_ySIZE() != 1 || rhs.get_zSIZE() != 1 ||
          rhs.get_bSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index x,y,z,b != 1 "
                        "cannot cast to ComplexField_AT");
      mem_alloc(rhs.get_aSIZE(), rhs.get_tSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_AT; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_AT; }
};

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in time (only tb index)
 * @brief complex double Field(t, b)
 * @brief t: time  coodination
 * @brief b: outer degree of freedom
 */
//--------------------------------------------------------------------------
class ComplexField_TB : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_TB   myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int t, const int b) {
      return m_field[aVtb(0,0,t,b)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int t, const int b) const {
      return m_field[aVtb(0,0,t,b)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_TB() {}
   ComplexField_TB(const myclass& other) : base(other) {}
   ComplexField_TB(const base&    other) : base(1, 1,
                                                other.get_tSIZE(),
                                                other.get_bSIZE()) {
      (*this) = other;
   }
   ComplexField_TB(const base& other, const int ia, const int ix, const int iy,
                   const int iz)
   : base(1, 1, other.get_tSIZE(), other.get_bSIZE()) {
      for (   int b=0; b<(*this).get_bSIZE(); b++)
         for (int t=0; t<(*this).get_tSIZE(); t++)
            (*this)(t,b) = other(ia,ix,iy,iz,t,b);
   }
   ComplexField_TB(const base& other, const int ia, const int ixyz)
   : base(1, 1, other.get_tSIZE(), other.get_bSIZE()) {
      for (   int b=0; b<(*this).get_bSIZE(); b++)
         for (int t=0; t<(*this).get_tSIZE(); t++)
            (*this)(t,b) = other(ia,ixyz,t,b);
   }
   ComplexField_TB(const int a_tSIZE, const int a_bSIZE)
   : base(1, 1, a_tSIZE, a_bSIZE) {}
   
   ~ComplexField_TB() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(1, 1, a_tSIZE, a_bSIZE);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_aSIZE() != 1 || rhs.get_xSIZE() != 1 || rhs.get_ySIZE() != 1 ||
          rhs.get_zSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index a,x,y,z != 1 "
                        "cannot cast to ComplexField_TB");
      mem_alloc(rhs.get_tSIZE(), rhs.get_bSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_TB; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_TB; }
};

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in 3-dimensional space (only xyz index)
 * @brief complex double Field(x, y, z)
 * @brief x,y,z: space coodination (x -> inner, z -> outer)
 */
//--------------------------------------------------------------------------
class ComplexField_XYZ : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_XYZ  myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int x, const int y, const int z) {
      return m_field[axyztb(0,x,y,z,0,0)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int x, const int y, const int z) const {
      return m_field[axyztb(0,x,y,z,0,0)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_XYZ() {}
   ComplexField_XYZ(const myclass& other) : base(other) {}
   ComplexField_XYZ(const base&    other) : base(1,
                                                 other.get_xSIZE(),
                                                 other.get_ySIZE(),
                                                 other.get_zSIZE(), 1, 1) {
      (*this) = other;
   }
   ComplexField_XYZ(const base& other, const int ia, const int it, const int ib)
   : base(1, other.get_xSIZE(), other.get_ySIZE(), other.get_zSIZE(), 1, 1) {
      for (int n=0; n<(*this).data_size(); n++) (*this)(n) = other(ia,n,it,ib);
   }
   ComplexField_XYZ(const int a_xSIZE, const int a_ySIZE, const int a_zSIZE)
   : base(1, a_xSIZE, a_ySIZE, a_zSIZE, 1, 1) {}
   ComplexField_XYZ(const int a_Lsize) : base(1, a_Lsize, 1, 1) {}
   
   ~ComplexField_XYZ() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_xSIZE, const int a_ySIZE, const int a_zSIZE) {
      base::mem_alloc(1, a_xSIZE, a_ySIZE, a_zSIZE, 1, 1);
   }
   void mem_alloc(const int a_Lsize) {
      base::mem_alloc(1, a_Lsize, 1, 1);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_aSIZE() != 1 || rhs.get_tSIZE() != 1 || rhs.get_bSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index a,t,b != 1 "
                        "cannot cast to ComplexField_XYZ");
      mem_alloc(rhs.get_xSIZE(), rhs.get_ySIZE(), rhs.get_zSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_XYZ; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_XYZ; }
};

//--------------------------------------------------------------------------
/**
 * @brief The spherical harmonics as a field in 3-dimension space
 */
//--------------------------------------------------------------------------
namespace sfunc {
   ComplexField_XYZ cfield_Ylm(const int, const int, const int);
}

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in 3-dimensional space
 * @brief complex double Field(a, x, y, z, b)
 * @brief a    : inner degree of freedom
 * @brief x,y,z: space coodination (x -> inner, z -> outer)
 * @brief b    : outer degree of freedom
 */
//--------------------------------------------------------------------------
class ComplexField_AXYZB : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE  base;
   typedef ComplexField_AXYZB myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int a, const int x, const int y, const int z,
                       const int b) {
      return m_field[axyztb(a,x,y,z,0,b)];
   }
   cdouble& operator()(const int a, const int xyz, const int b) {
      return m_field[aVtb(a,xyz,0,b)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int a, const int x, const int y, const int z,
                             const int b) const {
      return m_field[axyztb(a,x,y,z,0,b)];
   }
   const cdouble& operator()(const int a, const int xyz, const int b) const {
      return m_field[aVtb(a,xyz,0,b)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_AXYZB() {}
   ComplexField_AXYZB(const myclass& other) : base(other) {}
   ComplexField_AXYZB(const base&    other) : base(other.get_aSIZE(),
                                                   other.get_xSIZE(),
                                                   other.get_ySIZE(),
                                                   other.get_zSIZE(), 1,
                                                   other.get_bSIZE()) {
      (*this) = other;
   }
   ComplexField_AXYZB(const base& other, const int it)
   : base(other.get_aSIZE(), other.get_xSIZE(), other.get_ySIZE(),
          other.get_zSIZE(),                 1, other.get_bSIZE()) {
      for (      int b=0; b<(*this). get_bSIZE(); b++)
         for (   int n=0; n<(*this).data_Vsize(); n++)
            for (int a=0; a<(*this). get_aSIZE(); a++)
               (*this)(a,n,b) = other(a,n,it,b);
   }
   ComplexField_AXYZB(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                      const int a_zSIZE, const int a_bSIZE)
   : base(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, 1, a_bSIZE) {}
   ComplexField_AXYZB(const int a_aSIZE, const int a_Lsize, const int a_bSIZE)
   : base(a_aSIZE, a_Lsize, 1, a_bSIZE) {}
   
   ~ComplexField_AXYZB() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, 1, a_bSIZE);
   }
   void mem_alloc(const int a_aSIZE, const int a_Lsize, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_Lsize, 1, a_bSIZE);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_tSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index t != 1 "
                        "cannot cast to ComplexField_AXYZB");
      mem_alloc(rhs.get_aSIZE(), rhs.get_xSIZE(), rhs.get_ySIZE(),
                rhs.get_zSIZE(), rhs.get_bSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_AXYZB; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_AXYZB; }
};

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in 3-dimensional space (only axyz index)
 * @brief complex double Field(a, x, y, z)
 * @brief a    : inner degree of freedom
 * @brief x,y,z: space coodination (x -> inner, z -> outer)
 */
//--------------------------------------------------------------------------
class ComplexField_AXYZ : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_AXYZ myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int a, const int x, const int y, const int z) {
      return m_field[axyztb(a,x,y,z,0,0)];
   }
   cdouble& operator()(const int a, const int xyz) {
      return m_field[aVtb(a,xyz,0,0)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int a, const int x, const int y,
                             const int z) const {
      return m_field[axyztb(a,x,y,z,0,0)];
   }
   const cdouble& operator()(const int a, const int xyz) const {
      return m_field[aVtb(a,xyz,0,0)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_AXYZ() {}
   ComplexField_AXYZ(const myclass& other) : base(other) {}
   ComplexField_AXYZ(const base&    other) : base(other.get_aSIZE(),
                                                  other.get_xSIZE(),
                                                  other.get_ySIZE(),
                                                  other.get_zSIZE(), 1, 1) {
      (*this) = other;
   }
   ComplexField_AXYZ(const base& other, const int it, const int ib)
   : base(other.get_aSIZE(), other.get_xSIZE(), other.get_ySIZE(),
          other.get_zSIZE(),                 1,                 1) {
      for (   int n=0; n<(*this).data_Vsize(); n++)
         for (int a=0; a<(*this). get_aSIZE(); a++)
            (*this)(a,n) = other(a,n,it,ib);
   }
   ComplexField_AXYZ(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                     const int a_zSIZE)
   : base(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, 1, 1) {}
   ComplexField_AXYZ(const int a_aSIZE, const int a_Lsize)
   : base(a_aSIZE, a_Lsize, 1, 1) {}
   
   ~ComplexField_AXYZ() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, 1, 1);
   }
   void mem_alloc(const int a_aSIZE, const int a_Lsize) {
      base::mem_alloc(a_aSIZE, a_Lsize, 1, 1);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_tSIZE() != 1 || rhs.get_bSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index t,b != 1 "
                        "cannot cast to ComplexField_AXYZ");
      mem_alloc(rhs.get_aSIZE(), rhs.get_xSIZE(), rhs.get_ySIZE(),
                rhs.get_zSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_AXYZ; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_AXYZ; }
};

//--------------------------------------------------------------------------
/**
 * @brief The sub class of complex-field in 3-dimensional space (only xyzb index)
 * @brief complex double Field(x, y, z, b)
 * @brief x,y,z: space coodination (x -> inner, z -> outer)
 * @brief b    : outer degree of freedom
 */
//--------------------------------------------------------------------------
class ComplexField_XYZB : public ComplexField_BASE
{
public:
   typedef ComplexField_BASE base;
   typedef ComplexField_XYZB myclass;
   using   base::operator();
   using   base::operator =;
   using   base::operator+=;
   using   base::operator-=;
   using   base::operator*=;
   using   base::operator/=;
//============================== For writing =============================//
   cdouble& operator()(const int x, const int y, const int z, const int b) {
      return m_field[axyztb(0,x,y,z,0,b)];
   }
   cdouble& operator()(const int xyz, const int b) {
      return m_field[aVtb(0,xyz,0,b)];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int x, const int y, const int z,
                             const int b) const {
      return m_field[axyztb(0,x,y,z,0,b)];
   }
   const cdouble& operator()(const int xyz, const int b) const {
      return m_field[aVtb(0,xyz,0,b)];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_XYZB() {}
   ComplexField_XYZB(const myclass& other) : base(other) {}
   ComplexField_XYZB(const base&    other) : base(1,
                                                  other.get_xSIZE(),
                                                  other.get_ySIZE(),
                                                  other.get_zSIZE(), 1,
                                                  other.get_bSIZE()) {
      (*this) = other;
   }
   ComplexField_XYZB(const base& other, const int ia, const int it)
   : base(                1, other.get_xSIZE(), other.get_ySIZE(),
          other.get_zSIZE(),                 1, other.get_bSIZE()) {
      for (   int b=0; b<(*this). get_bSIZE(); b++)
         for (int n=0; n<(*this).data_Vsize(); n++)
            (*this)(n,b) = other(ia,n,it,b);
   }
   ComplexField_XYZB(const int a_xSIZE, const int a_ySIZE, const int a_zSIZE,
                     const int a_bSIZE)
   : base(1, a_xSIZE, a_ySIZE, a_zSIZE, 1, a_bSIZE) {}
   ComplexField_XYZB(const int a_Lsize, const int a_bSIZE)
   : base(1, a_Lsize, 1, a_bSIZE) {}
   
   ~ComplexField_XYZB() {}
//============================= For initialize ===========================//
   void mem_alloc(const int a_aSIZE, const int a_xSIZE, const int a_ySIZE,
                  const int a_zSIZE, const int a_tSIZE, const int a_bSIZE) {
      base::mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
   }
   void mem_alloc(const int a_xSIZE, const int a_ySIZE, const int a_zSIZE,
                  const int a_bSIZE) {
      base::mem_alloc(1, a_xSIZE, a_ySIZE, a_zSIZE, 1, a_bSIZE);
   }
   void mem_alloc(const int a_Lsize, const int a_bSIZE) {
      base::mem_alloc(1, a_Lsize, 1, a_bSIZE);
   }
//========================= Operator definitions ========================//
   myclass& operator =(const base &rhs) {
      if (rhs.get_tSIZE() != 1 || rhs.get_aSIZE() != 1)
         ERROR_COMMENTS("ComplexField_BASE with index a,t != 1 "
                        "cannot cast to ComplexField_XYZB");
      mem_alloc(rhs.get_xSIZE(), rhs.get_ySIZE(), rhs.get_zSIZE(),
                rhs.get_bSIZE());
      for (int i=0; i<(*this).data_size(); i++) (*this)(i) = rhs(i);
      return *this;
   }
//=========================== Several functions ==========================//
   int       data_type ()       { return DTYPE_COMPLEX_FIELD_XYZB; }
   const int data_type () const { return DTYPE_COMPLEX_FIELD_XYZB; }
};

template <> void STATISTICS<ComplexField_T> ::
output_data_err(const string, const double, const bool);

template <> void STATISTICS<ComplexField_XYZ> ::
output_data_err       (const string, const double, const bool);
template <> void STATISTICS<ComplexField_XYZ> ::
output_data_bin_reduce(const string, const double, const bool);

/*
template <> void STATISTICS<ComplexField_AXYZTB>:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_AXYZB >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_AXYZ  >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_XYZB  >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_XYZ   >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_ATB   >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_AT    >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_TB    >:: input_data_bin(const string);
template <> void STATISTICS<ComplexField_T     >:: input_data_bin(const string);

template <> void STATISTICS<ComplexField_AXYZTB>::output_data_bin(const string);
template <> void STATISTICS<ComplexField_AXYZB >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_AXYZ  >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_XYZB  >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_XYZ   >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_ATB   >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_AT    >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_TB    >::output_data_bin(const string);
template <> void STATISTICS<ComplexField_T     >::output_data_bin(const string);
*/

#endif
