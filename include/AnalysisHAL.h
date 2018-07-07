//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Definitions of the namespace for common variables & functions
 * @author  Takaya Miyamoto
 * @since   Mon Nov  7 16:14:42 JST 2016
 */
//--------------------------------------------------------------------------

#ifndef ANALYSIS_HAL_H
#define ANALYSIS_HAL_H

#include "CommonIncl.h"
#include "AnalysisHAL_const.h"

namespace anaHAL {
   //! For endian convert
   inline bool machine_is_little_endian() {
      int intONE = 1; return (*(char*)&intONE) ? true : false;
   }
   template <typename X>
   inline void endian_convert(X *DATA, const int DATA_size) {
      DEBUG_LOG
      const int dsize = sizeof(X);
      char     *tmp_c = new char[dsize];
      
      for (int k=0; k<DATA_size; k++) {
         for (int j=0; j<dsize; j++) tmp_c[j] = ((char*)&DATA[k])[j];
         for (int j=0; j<dsize; j++) ((char*)&DATA[k])[j] = tmp_c[dsize-j-1];
      }
      delete tmp_c;
   }
   //! Specialization of complex double
   template <>
   inline void   endian_convert(cdouble *DATA, const int DATA_size) {
      DEBUG_LOG; endian_convert((double*)DATA, 2 * DATA_size);
   }
   
   //! Convert the bool to string ( yes or no )
   inline string bool_to_str(const bool flg) {
      return flg ? "yes" : "no";
   }
   //! Convert the string ( yes or no ) to bool
   inline bool   str_to_bool(const string str) {
      return (str == "yes" || str == "y") ? true : false;
   }
   
   //! Convert index = idx[0] + iSIZE[0] * (idx[1] + iSIZE[1] * (idx[2] + ...
   //! to idx[0], idx[1], idx[2], ... (#.idx := Nidx)
   inline void convert_idx(const size_t  index,       int *idx,
                           const    int *iSIZE, const int Nidx) {
      size_t tmp_idx = index;
      for (int i=0; i<Nidx; i++) {
         idx[i]  =  tmp_idx         % iSIZE[i];
         tmp_idx = (tmp_idx-idx[i]) / iSIZE[i];
      } if (tmp_idx != 0) ERROR_COMMENTS("index convert failed.");
   }
   //! Convert index = 0,1,...,L to -L/2+1,...,L/2 for each element
   inline void convert_origin(const int *iorg, int *isft,
                              const int *iSIZE, const int Nidx) {
      for (int i=0; i<Nidx; i++)
         iorg[i] > iSIZE[i]/2 ? isft[i] = iorg[i]-iSIZE[i] : isft[i] = iorg[i];
   }
   
   //! Calculation of reduced #.data point
   int reduced_Ndata(const int, const int, const int);
   
   //! For statistics
   void make_mean_err(const  double*,  double&,  double&, const int, const bool);
   void make_mean_err(const cdouble*, cdouble&, cdouble&, const int, const bool);
   
   //! For rotation matrix
   void rot_matrix_init(int rot_matrix[24][4][4], int rot_character[5][24],
                        const int, const int, const int);
   inline void rot_matrix_init(int rot_matrix[24][4][4],
                               int rot_character[5][24], const int Lsize) {
      rot_matrix_init(rot_matrix, rot_character, Lsize, Lsize, Lsize);
   }
}

namespace anaHAL {
   //! The class for name list (gauge conf, hadron name, etc...)
   class NameList {
   protected:
      string *m_lists;
      int     m_Nlist;
      
   public:
      //================== For writing & reading ==================//
      string& operator()(const int index) {
         if (index >= m_Nlist) ERROR_COMMENTS("Index Overflow");
         return m_lists[index];
      }
      const string& operator()(const int index) const {
         if (index >= m_Nlist) ERROR_COMMENTS("Index Overflow");
         return m_lists[index];
      }
      //================= Constructor & Destructor =================//
      NameList();
      NameList(const NameList&);
      NameList(const string , const int rcolumn = 1);
      NameList(const string*, const int);
      ~NameList();
      //===================== For initialize =====================//
      void mem_alloc();
      void mem_alloc(const int);
      void mem_del  ();
      void set      (const string , const int rcolumn = 1);
      void set      (const string*, const int);
      //=================== Operator definitions ==================//
      NameList& operator=(const NameList &rhs) {
         mem_alloc(rhs.Nlist());
         for (int i=0; i<(*this).Nlist(); i++) (*this)(i) = rhs(i);
         return *this;
      }
      //================== Several function ==================//
      const int Nlist() const { return m_Nlist; }
   };
}

//! For special functions
namespace sfunc {
   cdouble Y_0_0 (const int, const int, const int);
   
   cdouble Y_1_m1(const int, const int, const int);
   cdouble Y_1_0 (const int, const int, const int);
   cdouble Y_1_p1(const int, const int, const int);
   
   cdouble Y_2_m2(const int, const int, const int);
   cdouble Y_2_m1(const int, const int, const int);
   cdouble Y_2_0 (const int, const int, const int);
   cdouble Y_2_p1(const int, const int, const int);
   cdouble Y_2_p2(const int, const int, const int);
   
   cdouble Y_3_m3(const int, const int, const int);
   cdouble Y_3_m2(const int, const int, const int);
   cdouble Y_3_m1(const int, const int, const int);
   cdouble Y_3_0 (const int, const int, const int);
   cdouble Y_3_p1(const int, const int, const int);
   cdouble Y_3_p2(const int, const int, const int);
   cdouble Y_3_p3(const int, const int, const int);
   
   cdouble Y_4_m4(const int, const int, const int);
   cdouble Y_4_m3(const int, const int, const int);
   cdouble Y_4_m2(const int, const int, const int);
   cdouble Y_4_m1(const int, const int, const int);
   cdouble Y_4_0 (const int, const int, const int);
   cdouble Y_4_p1(const int, const int, const int);
   cdouble Y_4_p2(const int, const int, const int);
   cdouble Y_4_p3(const int, const int, const int);
   cdouble Y_4_p4(const int, const int, const int);
}

#endif
