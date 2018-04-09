//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Common header file for statistics template class
 * @author  Takaya Miyamoto
 * @since   Mon Nov  7 16:14:42 JST 2016
 */
//--------------------------------------------------------------------------

#ifndef STATISTICS_TEMPLATE_H
#define STATISTICS_TEMPLATE_H

#include "AnalysisHAL.h"

//--------------------------------------------------------------------------
/**
 * @brief The template class for statistical data
 */
//--------------------------------------------------------------------------
template <class X> class STATISTICS
{
protected:
   X   *m_data;
   int  m_Ndata;
   
//============================ For inner index ===========================//
public:
//============================== For writing =============================//
   X& operator()(const int index) {
      return m_data[index];
   }
//============================== For reading =============================//
   const X& operator()(const int index) const {
      return m_data[index];
   }
//======================== Constructor & Destructor ======================//
   STATISTICS() {
      m_data  = NULL;
      m_Ndata = 0;
   }
   STATISTICS(const STATISTICS& other) {
      m_data  = NULL;
      m_Ndata = 0;
      (*this) = other;
   }
   STATISTICS(const int a_Ndata) {
      m_data  = NULL;
      m_Ndata = 0;
      mem_alloc(a_Ndata);
   }
   ~STATISTICS() {
      if (m_data != NULL) delete [] m_data;
   }
//============================= For initialize ===========================//
   void mem_alloc() {
      DEBUG_LOG
      if ((*this).Ndata() == 0)
         ERROR_COMMENTS("The data size has not been initialized.");
      if (m_data == NULL) m_data = new X[(*this).Ndata()];
   }
   void mem_alloc(const int a_Ndata) {
      DEBUG_LOG
      if ((*this).Ndata() != a_Ndata) {
         mem_del();
         m_Ndata = a_Ndata;
      }
      mem_alloc();
   }
   void mem_del() {
      DEBUG_LOG
      if (m_data != NULL) {
         delete [] m_data;
         m_data = NULL;
      }
   }
//============================ Operator define ===========================//
   STATISTICS& operator =(const STATISTICS &rhs) {
      mem_alloc(rhs.Ndata());
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n)  = rhs(n);
      return *this;
   }
   STATISTICS& operator+=(const STATISTICS &rhs) {
      if ((*this).Ndata() != rhs.Ndata()) ERROR_COMMENTS("Different data size.");
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) += rhs(n);
      return *this;
   }
   STATISTICS& operator-=(const STATISTICS &rhs) {
      if ((*this).Ndata() != rhs.Ndata()) ERROR_COMMENTS("Different data size.");
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) -= rhs(n);
      return *this;
   }
   STATISTICS& operator*=(const STATISTICS &rhs) {
      if ((*this).Ndata() != rhs.Ndata()) ERROR_COMMENTS("Different data size.");
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) *= rhs(n);
      return *this;
   }
   STATISTICS& operator/=(const STATISTICS &rhs) {
      if ((*this).Ndata() != rhs.Ndata()) ERROR_COMMENTS("Different data size.");
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) /= rhs(n);
      return *this;
   }
   STATISTICS& operator =(const X rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n)  = rhs;
      return *this;
   }
   STATISTICS& operator+=(const X rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) += rhs;
      return *this;
   }
   STATISTICS& operator-=(const X rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) -= rhs;
      return *this;
   }
   STATISTICS& operator*=(const X rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) *= rhs;
      return *this;
   }
   STATISTICS& operator/=(const X rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) /= rhs;
      return *this;
   }
   STATISTICS& operator =(const cdouble rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n)  = rhs;
      return *this;
   }
   STATISTICS& operator+=(const cdouble rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) += rhs;
      return *this;
   }
   STATISTICS& operator-=(const cdouble rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) -= rhs;
      return *this;
   }
   STATISTICS& operator*=(const cdouble rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) *= rhs;
      return *this;
   }
   STATISTICS& operator/=(const cdouble rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) /= rhs;
      return *this;
   }
   STATISTICS& operator =(const  double rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n)  = rhs;
      return *this;
   }
   STATISTICS& operator+=(const  double rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) += rhs;
      return *this;
   }
   STATISTICS& operator-=(const  double rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) -= rhs;
      return *this;
   }
   STATISTICS& operator*=(const  double rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) *= rhs;
      return *this;
   }
   STATISTICS& operator/=(const  double rhs) {
      for (int n=0; n<(*this).Ndata(); n++) (*this)(n) /= rhs;
      return *this;
   }
//============================ Operator helper ===========================//
   friend STATISTICS operator+(STATISTICS lhs, const STATISTICS &rhs) {
      return lhs += rhs;
   }
   friend STATISTICS operator-(STATISTICS lhs, const STATISTICS &rhs) {
      return lhs -= rhs;
   }
   friend STATISTICS operator*(STATISTICS lhs, const STATISTICS &rhs) {
      return lhs *= rhs;
   }
   friend STATISTICS operator/(STATISTICS lhs, const STATISTICS &rhs) {
      return lhs /= rhs;
   }
   friend STATISTICS operator+(STATISTICS lhs, const cdouble& rhs) {
      return lhs += rhs;
   }
   friend STATISTICS operator-(STATISTICS lhs, const cdouble& rhs) {
      return lhs -= rhs;
   }
   friend STATISTICS operator*(STATISTICS lhs, const cdouble& rhs) {
      return lhs *= rhs;
   }
   friend STATISTICS operator/(STATISTICS lhs, const cdouble& rhs) {
      return lhs /= rhs;
   }
   friend STATISTICS operator-(STATISTICS lhs, const  double& rhs) {
      return lhs -= rhs;
   }
   friend STATISTICS operator*(STATISTICS lhs, const  double& rhs) {
      return lhs *= rhs;
   }
   friend STATISTICS operator/(STATISTICS lhs, const  double& rhs) {
      return lhs /= rhs;
   }
   friend STATISTICS operator+(const cdouble& lhs, STATISTICS rhs) {
      return rhs += lhs;
   }
   friend STATISTICS operator-(const cdouble& lhs, STATISTICS rhs) {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend STATISTICS operator*(const cdouble& lhs, STATISTICS rhs) {
      return rhs *= lhs;
   }
   friend STATISTICS operator+(const  double& lhs, STATISTICS rhs) {
      return rhs += lhs;
   }
   friend STATISTICS operator-(const  double& lhs, STATISTICS rhs) {
      rhs *= (-1.0);
      return rhs += lhs;
   }
   friend STATISTICS operator*(const  double& lhs, STATISTICS rhs) {
      return rhs *= lhs;
   }
//=========================== Several functions ==========================//
   int       Ndata()       { return m_Ndata; }
   const int Ndata() const { return m_Ndata; }
   
   void make_JK_sample(const int);
   void make_mean_err ( double*,  double*, const bool);
   void make_mean_err (cdouble*, cdouble*, const bool);
   
   void  input_data_bin(const string);
   void output_data_bin(const string);
   
   void output_data_err       (const string, const double, const bool);
   void output_data_bin_reduce(const string, const double, const bool);
};

//--------------------------------------------------------------------------
/**
 * @brief The member function to make Jack-Knife samples
 */
//--------------------------------------------------------------------------
template <class X> void STATISTICS<X>::make_JK_sample(const int a_BinSIZE) {
   DEBUG_LOG
   if ((*this).Ndata()  % a_BinSIZE != 0)
      ERROR_COMMENTS("#.data is indivisible by bin-size.");
   if ((*this).Ndata() == a_BinSIZE)
      ERROR_COMMENTS("#.data == bin-size is not allowed.");
   
   int Nbin     = (*this).Ndata() / a_BinSIZE;
   X*  tmp_data = new X[Nbin];
   
   X sum((*this)(0));
   for (int i=1; i<(*this).Ndata(); i++) sum += (*this)(i);
   
   for (int i=0; i<Nbin; i++) {
      X sum_bin((*this)(a_BinSIZE*i));
      for (int b=1; b<a_BinSIZE; b++)
         sum_bin += (*this)(b + a_BinSIZE*i);
      tmp_data[i] = (sum - sum_bin) / double((*this).Ndata()-a_BinSIZE);
   }
   mem_alloc(Nbin);
   for (int i=0; i<Nbin; i++) (*this)(i) = tmp_data[i];
   delete [] tmp_data;
}

//--------------------------------------------------------------------------
/**
 * @brief The member function to make mean & error
 */
//--------------------------------------------------------------------------
template <class X>
void STATISTICS<X>::make_mean_err(double *mean, double *err,
                                  const bool is_jack_knife_data) {
   DEBUG_LOG
   double *tmp_data = new double[(*this).Ndata()];
   double  tmp_mean, tmp_err;
   
   for (int n=0; n<(*this)(0).data_size(); n++) {
      for (int i=0; i<(*this).Ndata(); i++) tmp_data[i] = (*this)(i)(n).real();
      anaHAL::make_mean_err(tmp_data, tmp_mean, tmp_err, (*this).Ndata(),
                            is_jack_knife_data);
      mean[n] = tmp_mean;
      err [n] = tmp_err;
   }
   delete [] tmp_data;
}
template <class X>
void STATISTICS<X>::make_mean_err(cdouble *mean, cdouble *err,
                                  const bool is_jack_knife_data) {
   DEBUG_LOG
   cdouble *tmp_data = new cdouble[(*this).Ndata()];
   cdouble  tmp_mean, tmp_err;
   
   for (int n=0; n<(*this)(0).data_size(); n++) {
      for (int i=0; i<(*this).Ndata(); i++) tmp_data[i] = (*this)(i)(n);
      anaHAL::make_mean_err(tmp_data, tmp_mean, tmp_err, (*this).Ndata(),
                            is_jack_knife_data);
      mean[n] = tmp_mean;
      err [n] = tmp_err;
   }
   delete [] tmp_data;
}

//--------------------------------------------------------------------------
/**
 * @brief Function for input/output binary data
 */
//--------------------------------------------------------------------------

namespace { int32_t MgNum_multi_complex_field = 94421393; }

//======================== miyamoto-format notation =======================//
//
//                        !! ALWAYS LITTLE ENDIAN !!
//
//      (1) Magic Number (94421393) (32-bit integer)
//      (2) aSIZE (inner  DoF)      (32-bit integer)
//      (3) xSIZE (x-coordination)  (32-bit integer)
//      (4) ySIZE (y-coordination)  (32-bit integer)
//      (5) zSIZE (z-coordination)  (32-bit integer)
//      (6) tSIZE (t-coordination)  (32-bit integer)
//      (7) bSIZE (outner DoF)      (32-bit integer)
//      (8) Ndata                   (32-bit integer)
//
//      (9) data ((2)*(3)*(4)*(5)*(6)*(7)*(8)* 2 * 64-bit double float)
//
//         -> data[re/im+2*(a+aS*(x+xS*(y+yS*(z+zS*(t+tS*(b+bS*n))))))]
//
//=========================================================================//

template <class X>
void STATISTICS<X>::input_data_bin(const string ifile_name) {
   DEBUG_LOG
   if ((*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AXYZTB &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AXYZB  &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AXYZ   &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_XYZB   &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_XYZ    &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_ATB    &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AT     &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_TB     &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_T)
      ERROR_COMMENTS("Invalid data type, cannot to output the data.");
   
   ifstream ifs(ifile_name.c_str(), ios::in | ios::binary);
   if (!ifs) ERROR_FOPEN(ifile_name.c_str());
   
   int32_t MgNum;
   ifs.read((char*)&MgNum, sizeof(int32_t));
   if (!anaHAL::machine_is_little_endian()) anaHAL::endian_convert(&MgNum, 1);
   if (MgNum != MgNum_multi_complex_field)
      ERROR_COMMENTS("This file may not be Miyamoto-format binary file.");
   
   int32_t tmp_a, tmp_x, tmp_y, tmp_z, tmp_t, tmp_b, tmp_n;
   ifs.read((char*)&tmp_a, sizeof(int32_t));
   ifs.read((char*)&tmp_x, sizeof(int32_t));
   ifs.read((char*)&tmp_y, sizeof(int32_t));
   ifs.read((char*)&tmp_z, sizeof(int32_t));
   ifs.read((char*)&tmp_t, sizeof(int32_t));
   ifs.read((char*)&tmp_b, sizeof(int32_t));
   ifs.read((char*)&tmp_n, sizeof(int32_t));
   
   if (!anaHAL::machine_is_little_endian()) {
      anaHAL::endian_convert(&tmp_a, 1);
      anaHAL::endian_convert(&tmp_x, 1);
      anaHAL::endian_convert(&tmp_y, 1);
      anaHAL::endian_convert(&tmp_z, 1);
      anaHAL::endian_convert(&tmp_t, 1);
      anaHAL::endian_convert(&tmp_b, 1);
      anaHAL::endian_convert(&tmp_n, 1);
   }
   mem_alloc(tmp_n);
   
   cdouble *tmp_data = new cdouble[tmp_a*tmp_x*tmp_y*tmp_z*tmp_t*tmp_b];
   
   for (int i=0; i<(*this).Ndata(); i++) {
      (*this)(i).mem_alloc(tmp_a, tmp_x, tmp_y, tmp_z, tmp_t, tmp_b);
      
      ifs.read((char*)&tmp_data[0], sizeof(cdouble) * (*this)(i).data_size());
      
      if (!anaHAL::machine_is_little_endian())
         anaHAL::endian_convert(tmp_data, (*this)(i).data_size());
      
      for (int n=0; n<(*this)(i).data_size(); n++) (*this)(i)(n) = tmp_data[n];
   }
   delete [] tmp_data;
   ifs.close();
}
template <class X>
void STATISTICS<X>::output_data_bin(const string ofile_name) {
   DEBUG_LOG
   if ((*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AXYZTB &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AXYZB  &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AXYZ   &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_XYZB   &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_XYZ    &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_ATB    &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_AT     &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_TB     &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_T)
      ERROR_COMMENTS("Invalid data type, cannot to output the data.");
   
   if ((*this).Ndata() == 0)
      ERROR_COMMENTS("The data size has not been initialized.");
   for (int i=0; i<(*this).Ndata(); i++)
      if ((*this)(i).data_size() == 0)
         ERROR_COMMENTS("The data size has not been initialized.");
   
   ofstream ofs(ofile_name.c_str(), ios::out | ios::binary);
   if (!ofs) ERROR_FOPEN(ofile_name.c_str());
   
   int32_t MgNum = MgNum_multi_complex_field;
   int32_t tmp_a = (*this)(0).get_aSIZE();
   int32_t tmp_x = (*this)(0).get_xSIZE();
   int32_t tmp_y = (*this)(0).get_ySIZE();
   int32_t tmp_z = (*this)(0).get_zSIZE();
   int32_t tmp_t = (*this)(0).get_tSIZE();
   int32_t tmp_b = (*this)(0).get_bSIZE();
   int32_t tmp_n = (*this).Ndata();
   
   if (!anaHAL::machine_is_little_endian()) {
      anaHAL::endian_convert( &MgNum, 1);
      anaHAL::endian_convert( &tmp_a, 1);
      anaHAL::endian_convert( &tmp_x, 1);
      anaHAL::endian_convert( &tmp_y, 1);
      anaHAL::endian_convert( &tmp_z, 1);
      anaHAL::endian_convert( &tmp_t, 1);
      anaHAL::endian_convert( &tmp_b, 1);
      anaHAL::endian_convert( &tmp_n, 1);
   }
   ofs.write((char*)&MgNum, sizeof(int32_t));
   ofs.write((char*)&tmp_a, sizeof(int32_t));
   ofs.write((char*)&tmp_x, sizeof(int32_t));
   ofs.write((char*)&tmp_y, sizeof(int32_t));
   ofs.write((char*)&tmp_z, sizeof(int32_t));
   ofs.write((char*)&tmp_t, sizeof(int32_t));
   ofs.write((char*)&tmp_b, sizeof(int32_t));
   ofs.write((char*)&tmp_n, sizeof(int32_t));
   
   for (int i=0; i<(*this).Ndata(); i++) {
      cdouble *tmp_data = new cdouble[(*this)(i).data_size()];
      
      for (int n=0; n<(*this)(i).data_size(); n++) tmp_data[n] = (*this)(i)(n);
      
      if (!anaHAL::machine_is_little_endian())
         anaHAL::endian_convert(tmp_data, (*this)(i).data_size());
      
      ofs.write((char*)&tmp_data[0], sizeof(cdouble) * (*this)(i).data_size());
      
      delete [] tmp_data;
   }
   ofs.close();
}

//--------------------------------------------------------------------------
/**
 * @brief Function to output data with error
 */
//--------------------------------------------------------------------------
template <class X>
void STATISTICS<X>::output_data_err(const string outfile_name,
                                    const double lattice_spacing,
                                    const bool   is_jack_knife_data) {
   DEBUG_LOG
   if ((*this)(0).data_type() != DTYPE_COMPLEX_FIELD_T &&
       (*this)(0).data_type() != DTYPE_COMPLEX_FIELD_XYZ)
      ERROR_COMMENTS("Invalid data type, cannot to output the data.");
}

//--------------------------------------------------------------------------
/**
 * @brief Function for binary data output with reduced data size
 */
//--------------------------------------------------------------------------
template <class X>
void STATISTICS<X>::output_data_bin_reduce(const string outfile_name,
                                           const double lattice_spacing,
                                           const bool   is_complex_data) {
   DEBUG_LOG
   if ((*this)(0).data_type() != DTYPE_COMPLEX_FIELD_XYZ)
      ERROR_COMMENTS("Invalid data type, cannot to output the data.");
}

#endif
