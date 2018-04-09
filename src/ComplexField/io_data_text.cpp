//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Implementations of input/output complex field for text data
 * @author  Takaya Miyamoto
 * @since   Mon Nov  7 16:14:42 JST 2016
 */
//--------------------------------------------------------------------------

#include <ComplexField_Base.h>

namespace  {
//! Read a text data [Nrow x Ncol]: return double(rdata[icol + Ncol * irow])
   inline void read_text(FILE* fp, double* rdata,
                         const int Nrow, const int Ncol) {
      char  buf[1024];
      char *c_ptr;
      for (int irow=0; irow<Nrow; irow++) {
         fgets(buf, sizeof(buf), fp);
         c_ptr = strtok(buf," \t\n\0");
         if (c_ptr == NULL)
            ERROR_COMMENTS("A blank line is not allowed for reading text data.");
         
         rdata[Ncol * irow] = atof(c_ptr);
         
         for (int icol=1; icol<Ncol; icol++) {
            c_ptr = strtok(NULL," \t\n\0");
            if (c_ptr == NULL)
               ERROR_COMMENTS("Unexpected #.column in text data");
            
            rdata[icol + Ncol * irow] = atof(c_ptr);
         }
      }
   }
//! Split file path to dirname and basename
   inline void split_fpath(const string &fpath, string &dname, string &fname) {
      char  tmp_s[1024]; strcpy(tmp_s, fpath.c_str());
      char *tmp_ptr;
      tmp_ptr = strrchr(tmp_s, '/'); *tmp_ptr = '\0';
      
      fname = tmp_ptr+1;
      dname = tmp_s;
   }
}

//--------------------------------------------------------------------------
/**
 * @brief Input/Output bare text file
 * @brief This function assumes that #.row = tSIZE and #.column = aSIZE * 2
 * @brief (* 2) means the real part and the imaginary part
 */
//--------------------------------------------------------------------------
void ComplexField_BASE::input_data_text(const string ifile_name,
                                        const int    idx_b) {
   DEBUG_LOG
   if ((*this).get_xSIZE() != 1 || (*this).get_ySIZE() != 1 ||
       (*this).get_zSIZE() != 1)
      ERROR_COMMENTS("Only index x,y,z = 1 is allowed for reading text data.");
   
   if ((*this).get_bSIZE() <= idx_b) ERROR_COMMENTS("Index b Overflow.");
   
   int l_aSIZE = (*this).get_aSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   cdouble *tmp_data = new cdouble[l_tSIZE * l_aSIZE];
   
   FILE* fp = fopen(ifile_name.c_str(), "r");
   if(fp == NULL) ERROR_FOPEN(ifile_name.c_str());
   
   read_text(fp, (double*)tmp_data, l_tSIZE, l_aSIZE * 2);
   fclose   (fp);
   
   for (int it=0; it<l_tSIZE; it++) for (int ia=0; ia<l_aSIZE; ia++)
      (*this)(ia,0,it,idx_b) = tmp_data[ia + l_aSIZE*it];
   
   delete [] tmp_data;
}

void ComplexField_BASE::output_data_text(const string ofile_name,
                                         const int    idx_b) {
   DEBUG_LOG
   if ((*this).get_xSIZE() != 1 || (*this).get_ySIZE() != 1 ||
       (*this).get_zSIZE() != 1)
      ERROR_COMMENTS("Only index x,y,z = 1 is allowed for writing text data.");
   
   if ((*this).get_bSIZE() <= idx_b) ERROR_COMMENTS("Index b Overflow.");
   
   int l_aSIZE = (*this).get_aSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   
   FILE* fp = fopen(ofile_name.c_str(), "w");
   if(fp == NULL) ERROR_FOPEN(ofile_name.c_str());
   
   for (int it=0; it<l_tSIZE; it++) {
      for (int ia=0; ia<l_aSIZE; ia++)
         fprintf(fp, "%1.16e %1.16e\t",
                 (*this)(ia,0,it,idx_b).real(),
                 (*this)(ia,0,it,idx_b).imag());
      fprintf(fp, "\n");
   }
   fclose(fp);
}

//--------------------------------------------------------------------------
/**
 * @brief Input/Output bare text file (For double variables)
 * @brief This function assumes that #.row = tSIZE and #.column = aSIZE
 */
//--------------------------------------------------------------------------
void ComplexField_BASE::input_data_text_real(const string ifile_name,
                                             const int    idx_b) {
   DEBUG_LOG
   if ((*this).get_xSIZE() != 1 || (*this).get_ySIZE() != 1 ||
       (*this).get_zSIZE() != 1)
      ERROR_COMMENTS("Only index x,y,z = 1 is allowed for reading text data.");
   
   if ((*this).get_bSIZE() <= idx_b) ERROR_COMMENTS("Index b Overflow.");
   
   int l_aSIZE = (*this).get_aSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   double *tmp_data = new double[l_tSIZE * l_aSIZE];
   
   FILE* fp = fopen(ifile_name.c_str(), "r");
   if(fp == NULL) ERROR_FOPEN(ifile_name.c_str());
   
   read_text(fp, tmp_data, l_tSIZE, l_aSIZE);
   fclose   (fp);
   
   for (int it=0; it<l_tSIZE; it++) for (int ia=0; ia<l_aSIZE; ia++)
      (*this)(ia,0,it,idx_b) = tmp_data[ia + l_aSIZE*it];
   
   delete [] tmp_data;
}

void ComplexField_BASE::output_data_text_real(const string ofile_name,
                                              const int    idx_b) {
   DEBUG_LOG
   if ((*this).get_xSIZE() != 1 || (*this).get_ySIZE() != 1 ||
       (*this).get_zSIZE() != 1)
      ERROR_COMMENTS("Only index x,y,z = 1 is allowed for writing text data.");
   
   if ((*this).get_bSIZE() <= idx_b) ERROR_COMMENTS("Index b Overflow.");
   
   int l_aSIZE = (*this).get_aSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   
   FILE* fp = fopen(ofile_name.c_str(), "w");
   if(fp == NULL) ERROR_FOPEN(ofile_name.c_str());
   
   for (int it=0; it<l_tSIZE; it++) {
      for (int ia=0; ia<l_aSIZE; ia++)
         fprintf(fp, "%1.16e ", (*this)(ia,0,it,idx_b).real());
      fprintf(fp, "\n");
   }
   fclose(fp);
}

//--------------------------------------------------------------------------
/**
 * @brief Input/Output text file for single hadron 2pt-correlator
 * @brief Automatic hadron-antihadron average can be done
 * @brief This function assumes that #.row = tSIZE and #.column = aSIZE * 2 + 1
 * @brief (+ 1) means the first(0th) column (info. of time-slice for corr-file)
 * @brief (* 2) means the real part and the imaginary part
 */
//--------------------------------------------------------------------------
void ComplexField_BASE::input_data_corr(const string ifile_name,
                                        const bool   fb_mean_flg,
                                        const int    idx_b) {
   DEBUG_LOG
   if ((*this).get_xSIZE() != 1 || (*this).get_ySIZE() != 1 ||
       (*this).get_zSIZE() != 1)
      ERROR_COMMENTS("Only index x,y,z = 1 is allowed for reading text data.");
   
   if ((*this).get_bSIZE() <= idx_b) ERROR_COMMENTS("Index b Overflow.");
   
   int l_aSIZE = (*this).get_aSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   double *tmp_data = new double[l_tSIZE * (l_aSIZE * 2 + 1)];
   
   FILE* fp = fopen(ifile_name.c_str(), "r");
   if(fp == NULL) ERROR_FOPEN(ifile_name.c_str());
   
   read_text(fp, tmp_data, l_tSIZE, (l_aSIZE * 2 + 1));
   fclose   (fp);
   
   for (int it=0; it<l_tSIZE; it++) {
      if (it != tmp_data[(l_aSIZE*2+1)*it])
         ERROR_COMMENTS("read a different time slice");
      
      for (int ia=0; ia<l_aSIZE; ia++)
      (*this)(ia,0,it,idx_b) = cdouble(tmp_data[1+0+2*ia+(l_aSIZE*2+1)*it],
                                       tmp_data[1+1+2*ia+(l_aSIZE*2+1)*it]);
   }
   DEBUG_COMMENTS("Correlator Readed.");
   
   if(fb_mean_flg) {   // Read anti-correlator files
      string dname, fname; split_fpath(ifile_name, dname, fname);
      fp = fopen((dname+"/anti"+fname).c_str(), "r");
      if(fp == NULL) ERROR_FOPEN((dname+"/anti"+fname).c_str());
      
      read_text(fp, tmp_data, l_tSIZE, (l_aSIZE * 2 + 1));
      fclose   (fp);
      
      for (int it=0; it<l_tSIZE; it++) {
         if (it != tmp_data[(l_aSIZE*2+1)*it])
            ERROR_COMMENTS("read a different time slice");
         
         for (int ia=0; ia<l_aSIZE; ia++)
            (*this)(ia,0,(l_tSIZE-it)%l_tSIZE,idx_b)
            += cdouble(tmp_data[1+0+2*ia+(l_aSIZE*2+1)*it],
                       tmp_data[1+1+2*ia+(l_aSIZE*2+1)*it]);
      }
      (*this) *= 0.5;
      DEBUG_COMMENTS("Anti-Correlator Readed.");
   }
   delete [] tmp_data;
}

void ComplexField_BASE::output_data_corr(const string ofile_name,
                                         const int    idx_b) {
   DEBUG_LOG
   if ((*this).get_xSIZE() != 1 || (*this).get_ySIZE() != 1 ||
       (*this).get_zSIZE() != 1)
      ERROR_COMMENTS("Only index x,y,z = 1 is allowed for writing text data.");
   
   if ((*this).get_bSIZE() <= idx_b) ERROR_COMMENTS("Index b Overflow.");
   
   int l_aSIZE = (*this).get_aSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   
   FILE* fp = fopen(ofile_name.c_str(), "w");
   if(fp == NULL) ERROR_FOPEN(ofile_name.c_str());
   
   for (int it=0; it<l_tSIZE; it++) {
      fprintf(fp, "%4d", it);
      for (int ia=0; ia<l_aSIZE; ia++)
         fprintf(fp, "\t%1.16e %1.16e",
                 (*this)(ia,0,it,idx_b).real(),
                 (*this)(ia,0,it,idx_b).imag());
      fprintf(fp, "\n");
   }
   fclose(fp);
}
