//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Definitions of the class for name list
 * @author  Takaya Miyamoto
 * @since   Mon Nov  7 16:14:42 JST 2016
 */
//--------------------------------------------------------------------------

#include <AnalysisHAL.h>

//================= Constructor & Destructor ==================//
anaHAL::NameList::NameList() {
   m_lists = NULL;
   m_Nlist = 0;
}
anaHAL::NameList::NameList(const anaHAL::NameList &other) {
   m_lists = NULL;
   m_Nlist = 0;
   (*this) = other;
}
anaHAL::NameList::NameList(const string ifile_name, const int rcolumn) {
   m_lists = NULL;
   m_Nlist = 0;
   set(ifile_name, rcolumn);
}
anaHAL::NameList::NameList(const string *a_lists, const int a_Nlist) {
   m_lists = NULL;
   m_Nlist = 0;
   set(a_lists, a_Nlist);
}
anaHAL::NameList::~NameList() {
   if (m_lists != NULL) delete [] m_lists;
}
//================== For initialize ==================//
void anaHAL::NameList::mem_alloc() {
   DEBUG_LOG
   if ((*this).Nlist() == 0)
      ERROR_COMMENTS("The data size has not been initialized.");
   if (m_lists == NULL) m_lists = new string[(*this).Nlist()];
}
void anaHAL::NameList::mem_alloc(const int a_Nlist) {
   DEBUG_LOG
   if ((*this).Nlist() != a_Nlist) {
      mem_del();
      m_Nlist = a_Nlist;
   }
   mem_alloc();
}
void anaHAL::NameList::mem_del() {
   DEBUG_LOG
   if (m_lists != NULL) {
      delete [] m_lists;
      m_lists = NULL;
   }
}
void anaHAL::NameList::set(const string ifile_name, const int rcolumn) {
   DEBUG_LOG
   char   buf[2048];
   char  *c_ptr;
   FILE *fp = fopen(ifile_name.c_str(), "r");
   if (fp == NULL) ERROR_FOPEN(ifile_name.c_str());
   
   int icount = 0;
   while (fgets(buf, sizeof(buf), fp) != NULL) {
      c_ptr = strtok(buf, " \t\n\0");
      if (c_ptr == NULL || c_ptr[0] == '#') continue;
      for (int i=1; i<rcolumn; i++) {
         if (c_ptr == NULL) continue;
         c_ptr = strtok(NULL, " \t\n\0");
      }
      if (c_ptr == NULL) continue;
      icount++;
   }
   if (icount == 0) ERROR_COMMENTS("The file is empty.");
   mem_alloc(icount);
   rewind   (fp);
   icount = 0;
   while (fgets(buf, sizeof(buf), fp) != NULL) {
      c_ptr = strtok(buf, " \t\n\0");
      if (c_ptr == NULL || c_ptr[0] == '#') continue;
      for (int i=1; i<rcolumn; i++) {
         if (c_ptr == NULL) continue;
         c_ptr = strtok(NULL, " \t\n\0");
      }
      if (c_ptr == NULL) continue;
      (*this)(icount) = c_ptr;
      icount++;
   }
}
void anaHAL::NameList::set(const string *a_lists, const int a_Nlist) {
   DEBUG_LOG
   mem_alloc(a_Nlist);
   for (int i=0; i<(*this).Nlist(); i++) (*this)(i) = a_lists[i];
}
