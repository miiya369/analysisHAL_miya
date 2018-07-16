//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Implementations of base class of complex-field
 * @author  Takaya Miyamoto
 * @since   Mon Oct 30 18:50:56 JST 2017
 */
//--------------------------------------------------------------------------

#include <ComplexField_Base.h>

//--------------------------------------------------------------------------
/**
 * @brief The base class of complex-field
 * @brief complex double Field(a, x, y, z, t, b)
 * @brief a      : inner degree of freedom
 * @brief x,y,z,t: space-time coodination (x -> inner, t -> outer)
 * @brief b      : outer degree of freedom
 */
//--------------------------------------------------------------------------

//======================== Constructor & Destructor ======================//
ComplexField_BASE::ComplexField_BASE() {
   m_field = NULL;
   m_aSIZE = 0;
   m_xSIZE = 0;
   m_ySIZE = 0;
   m_zSIZE = 0;
   m_tSIZE = 0;
   m_bSIZE = 0;
}
ComplexField_BASE::ComplexField_BASE(const ComplexField_BASE& other) {
   m_field = NULL;
   m_aSIZE = 0;
   m_xSIZE = 0;
   m_ySIZE = 0;
   m_zSIZE = 0;
   m_tSIZE = 0;
   m_bSIZE = 0;
   (*this) = other;
}
ComplexField_BASE::ComplexField_BASE(const int a_aSIZE, const int a_xSIZE,
                                     const int a_ySIZE, const int a_zSIZE,
                                     const int a_tSIZE, const int a_bSIZE) {
   m_field = NULL;
   m_aSIZE = 0;
   m_xSIZE = 0;
   m_ySIZE = 0;
   m_zSIZE = 0;
   m_tSIZE = 0;
   m_bSIZE = 0;
   mem_alloc(a_aSIZE, a_xSIZE, a_ySIZE, a_zSIZE, a_tSIZE, a_bSIZE);
}
ComplexField_BASE::ComplexField_BASE(const int a_aSIZE, const int a_Lsize,
                                     const int a_tSIZE, const int a_bSIZE) {
   m_field = NULL;
   m_aSIZE = 0;
   m_xSIZE = 0;
   m_ySIZE = 0;
   m_zSIZE = 0;
   m_tSIZE = 0;
   m_bSIZE = 0;
   mem_alloc(a_aSIZE, a_Lsize, a_tSIZE, a_bSIZE);
}
ComplexField_BASE::~ComplexField_BASE() {
   if (m_field != NULL) delete [] m_field;
}
//============================= For initialize ===========================//
void ComplexField_BASE::mem_alloc() {
   DEBUG_LOG
   if ((*this).data_size() == 0)
      ERROR_COMMENTS("The data size has not been initialized.");
   if (m_field == NULL) m_field = new cdouble[(*this).data_size()];
}
void ComplexField_BASE::mem_alloc(const int a_aSIZE, const int a_xSIZE,
                                  const int a_ySIZE, const int a_zSIZE,
                                  const int a_tSIZE, const int a_bSIZE) {
   DEBUG_LOG
   if ((*this).get_aSIZE() != a_aSIZE || (*this).get_xSIZE() != a_xSIZE ||
       (*this).get_ySIZE() != a_ySIZE || (*this).get_zSIZE() != a_zSIZE ||
       (*this).get_tSIZE() != a_tSIZE || (*this).get_bSIZE() != a_bSIZE) {
      mem_del();
      m_aSIZE = a_aSIZE;
      m_xSIZE = a_xSIZE;
      m_ySIZE = a_ySIZE;
      m_zSIZE = a_zSIZE;
      m_tSIZE = a_tSIZE;
      m_bSIZE = a_bSIZE;
   }
   mem_alloc();
}
void ComplexField_BASE::mem_alloc(const int a_aSIZE, const int a_Lsize,
                                  const int a_tSIZE, const int a_bSIZE) {
   DEBUG_LOG
   mem_alloc(a_aSIZE, a_Lsize, a_Lsize, a_Lsize, a_tSIZE, a_bSIZE);
}
void ComplexField_BASE::mem_del() {
   DEBUG_LOG
   if (m_field != NULL) {
      delete [] m_field;
      m_field = NULL;
   }
}
//============================ Operator definitions ===========================//
ComplexField_BASE& ComplexField_BASE::operator =(const ComplexField_BASE &rhs) {
   mem_alloc(rhs.get_aSIZE(), rhs.get_xSIZE(), rhs.get_ySIZE(),
             rhs.get_zSIZE(), rhs.get_tSIZE(), rhs.get_bSIZE());
   for (int n=0; n<(*this).data_size(); n++) (*this)(n)  = rhs(n);
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator+=(const ComplexField_BASE &rhs) {
   if ((*this).data_size() != rhs.data_size())
      ERROR_COMMENTS("Different data size.");
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) += rhs(n);
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator-=(const ComplexField_BASE &rhs) {
   if ((*this).data_size() != rhs.data_size())
      ERROR_COMMENTS("Different data size.");
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) -= rhs(n);
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator*=(const ComplexField_BASE &rhs) {
   if ((*this).data_size() != rhs.data_size())
      ERROR_COMMENTS("Different data size.");
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) *= rhs(n);
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator/=(const ComplexField_BASE &rhs) {
   if ((*this).data_size() != rhs.data_size())
      ERROR_COMMENTS("Different data size.");
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) /= rhs(n);
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator =(const cdouble &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n)  = rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator+=(const cdouble &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) += rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator-=(const cdouble &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) -= rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator*=(const cdouble &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) *= rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator/=(const cdouble &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) /= rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator =(const  double &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n)  = rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator+=(const  double &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) += rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator-=(const  double &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) -= rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator*=(const  double &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) *= rhs;
   return *this;
}
ComplexField_BASE& ComplexField_BASE::operator/=(const  double &rhs) {
   for (int n=0; n<(*this).data_size(); n++) (*this)(n) /= rhs;
   return *this;
}
//============================ Operator helper ===========================//
   
//=========================== Several functions ==========================//
void ComplexField_BASE::parity_average() {
   DEBUG_LOG
   int l_aSIZE = (*this).get_aSIZE();
   int l_xSIZE = (*this).get_xSIZE();
   int l_ySIZE = (*this).get_ySIZE();
   int l_zSIZE = (*this).get_zSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   int l_bSIZE = (*this).get_bSIZE();
   ComplexField_BASE tmp(l_aSIZE, l_xSIZE, l_ySIZE, l_zSIZE, l_tSIZE, l_bSIZE);
   for (               int b=0; b<l_bSIZE; b++)
      for (            int t=0; t<l_tSIZE; t++)
         for (         int z=0; z<l_zSIZE; z++)
            for (      int y=0; y<l_ySIZE; y++)
               for (   int x=0; x<l_xSIZE; x++)
                  for (int a=0; a<l_aSIZE; a++)
                     tmp(a,x,y,z,t,b) = (*this)(a,
                                                (l_xSIZE-x)%l_xSIZE,
                                                (l_ySIZE-y)%l_ySIZE,
                                                (l_zSIZE-z)%l_zSIZE,
                                                t,b);
   (*this) += tmp;
   (*this) *= 0.5;
}

cdouble ComplexField_BASE::average_space(const int a, const int t,
                                         const int b) const {
   DEBUG_LOG
   cdouble ret = 0.0;
   for (int n=0; n<(*this).data_Vsize(); n++) ret += (*this)(a,n,t,b);
   
   return (ret / double((*this).data_Vsize()));
}

cdouble ComplexField_BASE::lap(const int a, const int x, const int y,
                               const int z, const int t, const int b) const {
   int Xp = (x+1                    ) % (*this).get_xSIZE();
   int Xm = (x-1+(*this).get_xSIZE()) % (*this).get_xSIZE();
   int Yp = (y+1                    ) % (*this).get_ySIZE();
   int Ym = (y-1+(*this).get_ySIZE()) % (*this).get_ySIZE();
   int Zp = (z+1                    ) % (*this).get_zSIZE();
   int Zm = (z-1+(*this).get_zSIZE()) % (*this).get_zSIZE();
   return ((*this)(a, Xp,y ,z , t, b) + (*this)(a, Xm,y ,z , t, b) +
           (*this)(a, x ,Yp,z , t, b) + (*this)(a, x ,Ym,z , t, b) +
           (*this)(a, x ,y ,Zp, t, b) + (*this)(a, x ,y ,Zm, t, b) -
           (*this)(a, x ,y ,z , t, b) * 6.0);
}

ComplexField_BASE ComplexField_BASE::lap() const {
   DEBUG_LOG
   int l_aSIZE = (*this).get_aSIZE();
   int l_xSIZE = (*this).get_xSIZE();
   int l_ySIZE = (*this).get_ySIZE();
   int l_zSIZE = (*this).get_zSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   int l_bSIZE = (*this).get_bSIZE();
   ComplexField_BASE ret(l_aSIZE, l_xSIZE, l_ySIZE, l_zSIZE, l_tSIZE, l_bSIZE);
   for (               int b=0; b<l_bSIZE; b++)
      for (            int t=0; t<l_tSIZE; t++)
         for (         int z=0; z<l_zSIZE; z++)
            for (      int y=0; y<l_ySIZE; y++)
               for (   int x=0; x<l_xSIZE; x++)
                  for (int a=0; a<l_aSIZE; a++)
                     ret(a,x,y,z,t,b) = (*this).lap(a,x,y,z,t,b);
   return ret;
}

cdouble ComplexField_BASE::lap2(const int a, const int x, const int y,
				const int z, const int t, const int b) const {
   int Xp1 = (x+1                    ) % (*this).get_xSIZE();
   int Xm1 = (x-1+(*this).get_xSIZE()) % (*this).get_xSIZE();
   int Yp1 = (y+1                    ) % (*this).get_ySIZE();
   int Ym1 = (y-1+(*this).get_ySIZE()) % (*this).get_ySIZE();
   int Zp1 = (z+1                    ) % (*this).get_zSIZE();
   int Zm1 = (z-1+(*this).get_zSIZE()) % (*this).get_zSIZE();
   int Xp2 = (x+2                    ) % (*this).get_xSIZE();
   int Xm2 = (x-2+(*this).get_xSIZE()) % (*this).get_xSIZE();
   int Yp2 = (y+2                    ) % (*this).get_ySIZE();
   int Ym2 = (y-2+(*this).get_ySIZE()) % (*this).get_ySIZE();
   int Zp2 = (z+2                    ) % (*this).get_zSIZE();
   int Zm2 = (z-2+(*this).get_zSIZE()) % (*this).get_zSIZE();
   return (((*this)(a, Xp1,y  ,z  , t, b) + (*this)(a, Xm1,y  ,z  , t, b) +
            (*this)(a, x  ,Yp1,z  , t, b) + (*this)(a, x  ,Ym1,z  , t, b) +
            (*this)(a, x  ,y  ,Zp1, t, b) + (*this)(a, x  ,y  ,Zm1, t, b)) * 16.0
           -
           ((*this)(a, Xp2,y  ,z  , t, b) + (*this)(a, Xm2,y  ,z  , t, b) +
            (*this)(a, x  ,Yp2,z  , t, b) + (*this)(a, x  ,Ym2,z  , t, b) +
            (*this)(a, x  ,y  ,Zp2, t, b) + (*this)(a, x  ,y  ,Zm2, t, b))
           -
           (*this)(a, x ,y ,z , t, b) * 90.0) / 12.0;
}

ComplexField_BASE ComplexField_BASE::lap2() const {
   DEBUG_LOG
   int l_aSIZE = (*this).get_aSIZE();
   int l_xSIZE = (*this).get_xSIZE();
   int l_ySIZE = (*this).get_ySIZE();
   int l_zSIZE = (*this).get_zSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   int l_bSIZE = (*this).get_bSIZE();
   ComplexField_BASE ret(l_aSIZE, l_xSIZE, l_ySIZE, l_zSIZE, l_tSIZE, l_bSIZE);
   for (               int b=0; b<l_bSIZE; b++)
      for (            int t=0; t<l_tSIZE; t++)
         for (         int z=0; z<l_zSIZE; z++)
            for (      int y=0; y<l_ySIZE; y++)
               for (   int x=0; x<l_xSIZE; x++)
                  for (int a=0; a<l_aSIZE; a++)
                     ret(a,x,y,z,t,b) = (*this).lap2(a,x,y,z,t,b);
   return ret;
}

ComplexField_BASE ComplexField_BASE::rot_proj(int rep_type) const {
   DEBUG_LOG
   int l_aSIZE = (*this).get_aSIZE();
   int l_xSIZE = (*this).get_xSIZE();
   int l_ySIZE = (*this).get_ySIZE();
   int l_zSIZE = (*this).get_zSIZE();
   int l_tSIZE = (*this).get_tSIZE();
   int l_bSIZE = (*this).get_bSIZE();
   ComplexField_BASE ret(l_aSIZE, l_xSIZE, l_ySIZE, l_zSIZE, l_tSIZE, l_bSIZE);
   ComplexField_BASE tmp(      1, l_xSIZE, l_ySIZE, l_zSIZE,       1,       1);
   
   double factor = 0;
   if      (rep_type == ROT_REP_A1 || rep_type == ROT_REP_A2) factor = 1.0;
   else if (rep_type == ROT_REP_E)                            factor = 2.0;
   else if (rep_type == ROT_REP_T1 || rep_type == ROT_REP_T2) factor = 3.0;
   else ERROR_COMMENTS("Invalid representation type.");
   factor /= 24.0;
   
   int xyzSIZE[3] = {l_xSIZE, l_ySIZE, l_zSIZE};
   int rot_mat[24][4][4], rot_char[5][24], rot_xyz[3];
   anaHAL::rot_matrix_init(rot_mat, rot_char, l_xSIZE, l_ySIZE, l_zSIZE);
   
   for (      int b=0; b<l_bSIZE; b++)
      for (   int t=0; t<l_tSIZE; t++)
         for (int a=0; a<l_aSIZE; a++)
         {
            for (int n=0; n<tmp.data_size(); n++) tmp(n) = 0.0;
            
            for (int rot_type=0; rot_type < 24; rot_type++)
               for (         int z=0; z<l_zSIZE; z++)
                  for (      int y=0; y<l_ySIZE; y++)
                     for (   int x=0; x<l_xSIZE; x++) {
                        for (int k=0; k<3; k++)
                           rot_xyz[k] = (rot_mat[rot_type][k][0] * x +
                                         rot_mat[rot_type][k][1] * y +
                                         rot_mat[rot_type][k][2] * z +
                                         rot_mat[rot_type][k][3] * 1
                                         ) % xyzSIZE[k];
                        tmp(0,rot_xyz[0],rot_xyz[1],rot_xyz[2],0,0) +=
                        ((*this)(a,x,y,z,t,b) *
                         double(rot_char[rep_type][rot_type]));
                     }
            tmp *= factor;
            for (int n=0; n<tmp.data_size(); n++) ret(a,n,t,b) = tmp(n);
         }
   return ret;
}
