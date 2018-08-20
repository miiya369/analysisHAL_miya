//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Potential
 * @brief   Implementations for coupled-channel effective central potentials
 * @author  Takaya Miyamoto
 * @since   Tue Feb  2 18:57:38 JST 2016
 */
//--------------------------------------------------------------------------

#include <Potential.h>
#include <ComplexMatrix.h>

#define ij (j+a_Ndim*i)

void Potential::calc_CCpotential_T1(      ComplexField_BASE *V,
                                    const ComplexField_BASE *Rm,
                                    const ComplexField_BASE *R0,
                                    const  double *m1,  const double *m2,
                                    const cdouble *fac, const int a_Ndim) {
   DEBUG_LOG
   ComplexField_BASE *KR = new ComplexField_BASE[a_Ndim*a_Ndim];
   
   for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
      V[ij].mem_alloc(R0[ij].get_aSIZE(),R0[ij].get_xSIZE(),R0[ij].get_ySIZE(),
                      R0[ij].get_zSIZE(),R0[ij].get_tSIZE(),R0[ij].get_bSIZE());
      
      KR[ij] = Potential::get_potential_T1_nume(Rm[ij], R0[ij], m1[i], m2[i]);
   }
   cmatrix Rmat(a_Ndim);
   cmatrix Kmat(a_Ndim);
   cmatrix Vmat(a_Ndim);
   
   for (int n=0; n<V[0].data_size(); n++) {
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
         Rmat(i,j) = R0[ij](n);
         Kmat(i,j) = KR[ij](n);
      }
      Vmat = Kmat * Rmat.inverce();
      
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++)
         V[ij](n) = Vmat(i,j) / fac[ij];
   }
   delete [] KR;
}
void Potential::calc_CCpotential_T1(      ComplexField_BASE *V,
                                    const ComplexField_BASE *Rm,
                                    const ComplexField_BASE *R0,
                                    const ComplexField_BASE *Rp,
                                    const  double *m1,  const double *m2,
                                    const cdouble *fac, const int a_Ndim) {
   DEBUG_LOG
   ComplexField_BASE *KR = new ComplexField_BASE[a_Ndim*a_Ndim];
   
   for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
      V[ij].mem_alloc(R0[ij].get_aSIZE(),R0[ij].get_xSIZE(),R0[ij].get_ySIZE(),
                      R0[ij].get_zSIZE(),R0[ij].get_tSIZE(),R0[ij].get_bSIZE());
      
      KR[ij] = Potential::get_potential_T1_nume(Rm[ij], R0[ij], Rp[ij],
                                                m1[i] , m2[i]);
   }
   cmatrix Rmat(a_Ndim);
   cmatrix Kmat(a_Ndim);
   cmatrix Vmat(a_Ndim);
   
   for (int n=0; n<V[0].data_size(); n++) {
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
         Rmat(i,j) = R0[ij](n);
         Kmat(i,j) = KR[ij](n);
      }
      Vmat = Kmat * Rmat.inverce();
      
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++)
         V[ij](n) = Vmat(i,j) / fac[ij];
   }
   delete [] KR;
}
void Potential::calc_CCpotential_T2(      ComplexField_BASE *V,
                                    const ComplexField_BASE *Rm,
                                    const ComplexField_BASE *R0,
                                    const ComplexField_BASE *Rp,
                                    const  double *m1,  const double *m2,
                                    const cdouble *fac, const int a_Ndim) {
   DEBUG_LOG
   ComplexField_BASE *KR = new ComplexField_BASE[a_Ndim*a_Ndim];
   
   for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
      V[ij].mem_alloc(R0[ij].get_aSIZE(),R0[ij].get_xSIZE(),R0[ij].get_ySIZE(),
                      R0[ij].get_zSIZE(),R0[ij].get_tSIZE(),R0[ij].get_bSIZE());
      
      KR[ij] = Potential::get_potential_T2_nume(Rm[ij], R0[ij], Rp[ij],
                                                m1[i] , m2[i]);
   }
   cmatrix Rmat(a_Ndim);
   cmatrix Kmat(a_Ndim);
   cmatrix Vmat(a_Ndim);
   
   for (int n=0; n<V[0].data_size(); n++) {
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
         Rmat(i,j) = R0[ij](n);
         Kmat(i,j) = KR[ij](n);
      }
      Vmat = Kmat * Rmat.inverce();
      
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++)
         V[ij](n) = Vmat(i,j) / fac[ij];
   }
   delete [] KR;
}

#undef ij
