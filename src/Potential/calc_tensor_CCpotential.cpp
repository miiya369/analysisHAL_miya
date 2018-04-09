//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Potential
 * @brief   Implementations for coupled-channel central & tensor potentials
 * @author  Takaya Miyamoto
 * @since   Tue Feb  2 18:57:38 JST 2016
 */
//--------------------------------------------------------------------------

#include <Potential.h>
#include <ComplexMatrix.h>

#define ij (j+a_Ndim*i)

void Potential::calc_tensor_CCpot_T2(      ComplexField_BASE *V_C,
                                           ComplexField_BASE *V_T,
                                     const ComplexField_BASE *Rm,
                                     const ComplexField_BASE *R0,
                                     const ComplexField_BASE *Rp,
                                     const  double *m1,  const double *m2,
                                     const cdouble *fac, const int a_Ndim) {
   DEBUG_LOG
   int l_xSize = R0[0].get_xSIZE();
   int l_ySize = R0[0].get_ySIZE();
   int l_zSize = R0[0].get_zSIZE();
   int l_tSize = R0[0].get_tSIZE();
   for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
      V_C[ij].mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
      V_T[ij].mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
      V_C[ij] = 0.0; V_T[ij] = 0.0;
   }
   
   ComplexField_BASE *R0_S = new ComplexField_BASE[a_Ndim*a_Ndim];
   ComplexField_BASE *R0_D = new ComplexField_BASE[a_Ndim*a_Ndim];
   ComplexField_BASE *KR_S = new ComplexField_BASE[a_Ndim*a_Ndim];
   ComplexField_BASE *KR_D = new ComplexField_BASE[a_Ndim*a_Ndim];
   ComplexField_BASE *SR_S = new ComplexField_BASE[a_Ndim*a_Ndim];
   ComplexField_BASE *SR_D = new ComplexField_BASE[a_Ndim*a_Ndim];
   ComplexField_BASE  Rm_tmp, R0_tmp, Rp_tmp, KR, SR;
   
   cmatrix R0Smat(a_Ndim), R0Dmat(a_Ndim);
   cmatrix KRSmat(a_Ndim), KRDmat(a_Ndim);
   cmatrix SRSmat(a_Ndim), SRDmat(a_Ndim);
   cmatrix  VCmat(a_Ndim),  VTmat(a_Ndim);
   
   for (int src_spin1_z = -1; src_spin1_z <= +1; src_spin1_z++) {
      for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
         Rm_tmp = Rm[ij].src_spin_proj(HH_OctOct, 1, src_spin1_z);
         R0_tmp = R0[ij].src_spin_proj(HH_OctOct, 1, src_spin1_z);
         Rp_tmp = Rp[ij].src_spin_proj(HH_OctOct, 1, src_spin1_z);
         
         R0_S[ij] = R0_tmp.rot_proj(ROT_REP_A1);
         R0_D[ij] = Potential::get_YDstar_Psi(R0_tmp-R0_S[ij], src_spin1_z);
         
         KR       = Potential::get_potential_T2_nume(Rm_tmp, R0_tmp,
                                                     Rp_tmp, m1[i] , m2[i]);
         KR_S[ij] = KR.rot_proj(ROT_REP_A1);
         KR_D[ij] = Potential::get_YDstar_Psi(KR-KR_S[ij], src_spin1_z);
         
         SR       = Potential::get_S12_Psi(R0_tmp);
         SR_S[ij] = SR.rot_proj(ROT_REP_A1);
         SR_D[ij] = Potential::get_YDstar_Psi(SR-SR_S[ij], src_spin1_z);
         
         R0_S[ij] = R0_S[ij].snk_spin_proj(HH_OctOct, 1, src_spin1_z);
         KR_S[ij] = KR_S[ij].snk_spin_proj(HH_OctOct, 1, src_spin1_z);
         SR_S[ij] = SR_S[ij].snk_spin_proj(HH_OctOct, 1, src_spin1_z);
         
         R0_D[ij] = (R0_D[ij].snk_spin_proj(HH_OctOct, 1,+1) +
                     R0_D[ij].snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
         KR_D[ij] = (KR_D[ij].snk_spin_proj(HH_OctOct, 1,+1) +
                     KR_D[ij].snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
         SR_D[ij] = (SR_D[ij].snk_spin_proj(HH_OctOct, 1,+1) +
                     SR_D[ij].snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      }
      for (int n=0; n<V_C[0].data_size(); n++) {
         for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
            R0Smat(i,j) = R0_S[ij](n);
            R0Dmat(i,j) = R0_D[ij](n);
            KRSmat(i,j) = KR_S[ij](n);
            KRDmat(i,j) = KR_D[ij](n);
            SRSmat(i,j) = SR_S[ij](n);
            SRDmat(i,j) = SR_D[ij](n);
         }
         cmatrix R0Smat_inv = R0Smat.inverce();
         cmatrix R0Dmat_inv = R0Dmat.inverce();
         cmatrix SRSmat_inv = SRSmat.inverce();
         cmatrix SRDmat_inv = SRDmat.inverce();
         
         VCmat = ((KRSmat*SRSmat_inv - KRDmat*SRDmat_inv) *
                  (R0Smat*SRSmat_inv - R0Dmat*SRDmat_inv).inverce());
         VTmat = ((KRSmat*R0Smat_inv - KRDmat*R0Dmat_inv) *
                  (SRSmat*R0Smat_inv - SRDmat*R0Dmat_inv).inverce());
         
         for (int i=0; i<a_Ndim; i++) for (int j=0; j<a_Ndim; j++) {
            V_C[ij](n) += VCmat(i,j) / (3.0 * fac[ij]);
            V_T[ij](n) += VTmat(i,j) / (3.0 * fac[ij]);
         }
      }
   }
   delete [] R0_S;
   delete [] R0_D;
   delete [] KR_S;
   delete [] KR_D;
   delete [] SR_S;
   delete [] SR_D;
}

#undef ij
