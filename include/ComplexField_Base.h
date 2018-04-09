//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Definitions for base class of complex-field
 * @author  Takaya Miyamoto
 * @since   Mon Oct 30 18:50:56 JST 2017
 */
//--------------------------------------------------------------------------

#ifndef COMPLEX_FIELD_BASE_H
#define COMPLEX_FIELD_BASE_H

#include "AnalysisHAL.h"

//--------------------------------------------------------------------------
/**
 * @brief The base class of complex-field
 * @brief complex double Field(a, x, y, z, t, b)
 * @brief a      : inner degree of freedom
 * @brief x,y,z,t: space-time coodination (x -> inner, t -> outer)
 * @brief b      : outer degree of freedom
 */
//--------------------------------------------------------------------------
class ComplexField_BASE
{
protected:
   cdouble *m_field;
   int      m_xSIZE; // size of x-coordination
   int      m_ySIZE; // size of y-coordination
   int      m_zSIZE; // size of z-coordination
   int      m_tSIZE; // size of t-coordination
   int      m_aSIZE; // size of inner degree of freedom
   int      m_bSIZE; // size of outer degree of freedom
   
//============================ For inner index ===========================//
   int axyztb(const int a,
              const int x, const int y, const int z, const int t,
              const int b) const {
      return a + m_aSIZE *(x + m_xSIZE *( y + m_ySIZE *
                                         (z + m_zSIZE *(t + m_tSIZE * b))));
   }
   int aVtb(const int a, const int xyz, const int t, const int b) const {
      return a + m_aSIZE *(xyz + m_xSIZE*m_ySIZE*m_zSIZE *(t + m_tSIZE * b));
   }
public:
   typedef ComplexField_BASE myclass;
//============================== For writing =============================//
   cdouble& operator()(const int a,
                       const int x, const int y, const int z, const int t,
                       const int b) {
      return m_field[axyztb(a,x,y,z,t,b)];
   }
   cdouble& operator()(const int a, const int xyz, const int t, const int b) {
      return m_field[aVtb(a,xyz,t,b)];
   }
   cdouble& operator()(const int index) {
      return m_field[index];
   }
//============================== For reading =============================//
   const cdouble& operator()(const int a,
                             const int x, const int y, const int z, const int t,
                             const int b) const {
      return m_field[axyztb(a,x,y,z,t,b)];
   }
   const cdouble& operator()(const int a, const int xyz,
                             const int t, const int b) const {
      return m_field[aVtb(a,xyz,t,b)];
   }
   const cdouble& operator()(const int index) const {
      return m_field[index];
   }
//======================== Constructor & Destructor ======================//
   ComplexField_BASE();
   ComplexField_BASE(const myclass&);
   ComplexField_BASE(const int, const int, const int,
                     const int, const int, const int);
   ComplexField_BASE(const int, const int, const int, const int);
   ~ComplexField_BASE();
//============================= For initialize ===========================//
   void mem_alloc();
   void mem_alloc(const int, const int, const int,
                  const int, const int, const int);
   void mem_alloc(const int, const int, const int, const int);
   void mem_del  ();
//========================= Operator definitions ========================//
   myclass& operator =(const myclass&);
   myclass& operator+=(const myclass&);
   myclass& operator-=(const myclass&);
   myclass& operator*=(const myclass&);
   myclass& operator/=(const myclass&);
   myclass& operator =(const cdouble&);
   myclass& operator+=(const cdouble&);
   myclass& operator-=(const cdouble&);
   myclass& operator*=(const cdouble&);
   myclass& operator/=(const cdouble&);
   myclass& operator =(const  double&);
   myclass& operator+=(const  double&);
   myclass& operator-=(const  double&);
   myclass& operator*=(const  double&);
   myclass& operator/=(const  double&);
//============================ Operator helper ===========================//
   friend myclass operator+(myclass lhs, const myclass &rhs) {
      return lhs += rhs;
   }
   friend myclass operator-(myclass lhs, const myclass &rhs) {
      return lhs -= rhs;
   }
   friend myclass operator*(myclass lhs, const myclass &rhs) {
      return lhs *= rhs;
   }
   friend myclass operator/(myclass lhs, const myclass &rhs) {
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
   friend myclass operator/(const cdouble &lhs, myclass rhs) {
      return rhs.inverce() *= lhs;
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
   friend myclass operator/(const  double &lhs, myclass rhs) {
      return rhs.inverce() *= lhs;
   }
//=========================== Several functions ==========================//
   int data_size  () { return m_aSIZE*m_xSIZE*m_ySIZE*m_zSIZE*m_tSIZE*m_bSIZE; }
   int data_VTsize() { return         m_xSIZE*m_ySIZE*m_zSIZE*m_tSIZE        ; }
   int data_Vsize () { return         m_xSIZE*m_ySIZE*m_zSIZE                ; }
   int data_type  () { return DTYPE_COMPLEX_FIELD_AXYZTB; }
   int get_xSIZE  () { return m_xSIZE; }
   int get_ySIZE  () { return m_ySIZE; }
   int get_zSIZE  () { return m_zSIZE; }
   int get_tSIZE  () { return m_tSIZE; }
   int get_aSIZE  () { return m_aSIZE; }
   int get_bSIZE  () { return m_bSIZE; }
   
   const int data_size  () const
   { return m_aSIZE*m_xSIZE*m_ySIZE*m_zSIZE*m_tSIZE*m_bSIZE; }
   const int data_VTsize() const
   { return         m_xSIZE*m_ySIZE*m_zSIZE*m_tSIZE        ; }
   const int data_Vsize () const
   { return         m_xSIZE*m_ySIZE*m_zSIZE                ; }
   const int data_type  () const { return DTYPE_COMPLEX_FIELD_AXYZTB; }
   const int get_xSIZE  () const { return m_xSIZE; }
   const int get_ySIZE  () const { return m_ySIZE; }
   const int get_zSIZE  () const { return m_zSIZE; }
   const int get_tSIZE  () const { return m_tSIZE; }
   const int get_aSIZE  () const { return m_aSIZE; }
   const int get_bSIZE  () const { return m_bSIZE; }
   
   void  input_data_miya(const string);
   void  input_data_bin (const string);
   void  input_data_comp(const string, const int it    = 0);
   void  input_data_text(const string, const int idx_b = 0);
   void  input_data_text_real(const string, const int idx_b = 0);
   void  input_data_corr(const string, const bool, const int idx_b = 0);
   
   void output_data_miya(const string);
   void output_data_bin (const string);
   void output_data_comp(const string, const string, const int it = 0);
   void output_data_text(const string, const int idx_b = 0);
   void output_data_text_real(const string, const int idx_b = 0);
   void output_data_corr(const string, const int idx_b = 0);
   
   myclass inverce() {
      DEBUG_LOG
      myclass ret((*this).get_aSIZE(), (*this).get_xSIZE(),
                  (*this).get_ySIZE(), (*this).get_zSIZE(),
                  (*this).get_tSIZE(), (*this).get_bSIZE());
      for (int n=0; n<(*this).data_size(); n++) ret(n) = 1.0 / (*this)(n);
      return ret;
   }
   myclass conj() {
      DEBUG_LOG
      myclass ret((*this).get_aSIZE(), (*this).get_xSIZE(),
                  (*this).get_ySIZE(), (*this).get_zSIZE(),
                  (*this).get_tSIZE(), (*this).get_bSIZE());
      for (int n=0; n<(*this).data_size(); n++) ret(n) = std::conj((*this)(n));
      return ret;
   }
   
   void    parity_average();
   cdouble average_space (const int, const int, const int) const;
   
   cdouble lap(const int, const int, const int,
               const int, const int, const int) const;
   myclass lap() const;
   
   myclass      rot_proj(const int) const;
   myclass src_spin_proj(const int, const int, const int) const;
   myclass snk_spin_proj(const int, const int, const int) const;
   myclass     spin_proj(const int, const int, const int,
                         const int, const int, const int) const;
};

#endif
