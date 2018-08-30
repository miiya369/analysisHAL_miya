//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Potential
 * @brief   Implementations for single-channel central & tensor potentials
 * @author  Takaya Miyamoto
 * @since   Tue Feb  2 18:57:38 JST 2016
 */
//--------------------------------------------------------------------------

#include <Potential.h>

//! Return [S_{12} Psi]_{alpha,beta} for tensor force (Octet baryon)
//! Require that the alpha,beta(#.spin index of sink) == 2 ((#.inner index == 4))
ComplexField_BASE Potential::get_S12_Psi(const ComplexField_BASE& Psi) {
   DEBUG_LOG
   int l_aSize = Psi.get_aSIZE();
   int l_xSize = Psi.get_xSIZE();
   int l_ySize = Psi.get_ySIZE();
   int l_zSize = Psi.get_zSIZE();
   int l_tSize = Psi.get_tSIZE();
   int l_bSize = Psi.get_bSIZE();
   
   if (l_aSize != 4) ERROR_COMMENTS("Only a=4 is allowed for S12_Psi.");
   ComplexField_BASE Ret(l_aSize, l_xSize, l_ySize, l_zSize, l_tSize, l_bSize);
   
   int XYZ[3], XYZsize[3] = {l_xSize, l_ySize, l_zSize};
   for (      int x=0; x<l_xSize; x++)
      for (   int y=0; y<l_ySize; y++)
         for (int z=0; z<l_zSize; z++) {
            int xyz[3] = {x,y,z};
            anaHAL::convert_origin(xyz, XYZ, XYZsize, 3);
            int X = XYZ[0], Y = XYZ[1], Z = XYZ[2];
            
            for (   int t=0; t<l_tSize; t++)
               for (int b=0; b<l_bSize; b++) {
                  Ret(0,x,y,z,t,b) =
                  (
                   sqrt(2.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(0,x,y,z,t,b)  +
                   sqrt(3.0) * sfunc::Y_2_m1(X,Y,Z) *(Psi(1,x,y,z,t,b)  +
                                                      Psi(2,x,y,z,t,b)) +
                   sqrt(12.0)* sfunc::Y_2_m2(X,Y,Z) * Psi(3,x,y,z,t,b)
                   ) * sqrt(8.0*PI/5.0);
                  
                  Ret(1,x,y,z,t,b) =
                  (
                   sqrt(3.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(0,x,y,z,t,b)  +
                   sqrt(2.0) * sfunc::Y_2_0 (X,Y,Z) *(Psi(1,x,y,z,t,b)  +
                                                      Psi(2,x,y,z,t,b)) +
                   sqrt(3.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(3,x,y,z,t,b)
                   ) * (-sqrt(8.0*PI/5.0));
            
                  Ret(2,x,y,z,t,b) = Ret(1,x,y,z,t,b);
                  
                  Ret(3,x,y,z,t,b) =
                  (
                   sqrt(12.0)* sfunc::Y_2_p2(X,Y,Z) * Psi(0,x,y,z,t,b)  +
                   sqrt(3.0) * sfunc::Y_2_p1(X,Y,Z) *(Psi(1,x,y,z,t,b)  +
                                                      Psi(2,x,y,z,t,b)) +
                   sqrt(2.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(3,x,y,z,t,b)
                   ) * sqrt(8.0*PI/5.0);
               }
         }
   return Ret;
}
//! Return [YD^star * Psi]_{alpha,beta} for tensor force (Octet baryon)
//! Require that the alpha,beta(#.spin index of sink) == 2 ((#.inner index == 4))
ComplexField_BASE Potential::get_YDstar_Psi(const ComplexField_BASE& Psi,
                                            const int src_spin1_z) {
   DEBUG_LOG
   int l_aSize = Psi.get_aSIZE();
   int l_xSize = Psi.get_xSIZE();
   int l_ySize = Psi.get_ySIZE();
   int l_zSize = Psi.get_zSIZE();
   int l_tSize = Psi.get_tSIZE();
   int l_bSize = Psi.get_bSIZE();
   
   if (l_aSize != 4) ERROR_COMMENTS("Only a=4 is allowed for S12_Psi.");
   ComplexField_BASE Ret(l_aSize, l_xSize, l_ySize, l_zSize, l_tSize, l_bSize);
   
   int XYZ[3], XYZsize[3] = {l_xSize, l_ySize, l_zSize};
   if (src_spin1_z == -1)
      for (      int x=0; x<l_xSize; x++)
         for (   int y=0; y<l_ySize; y++)
            for (int z=0; z<l_zSize; z++) {
               int xyz[3] = {x,y,z};
               anaHAL::convert_origin(xyz, XYZ, XYZsize, 3);
               int X = XYZ[0], Y = XYZ[1], Z = XYZ[2];
               
               for (   int t=0; t<l_tSize; t++)
                  for (int b=0; b<l_bSize; b++) {
                     Ret(0,x,y,z,t,b) =
                     sqrt(1.0/10.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(0,x,y,z,t,b) +
                     sqrt(3.0/20.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(2,x,y,z,t,b);
                     
                     Ret(1,x,y,z,t,b) =
                     sqrt(1.0/10.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(1,x,y,z,t,b) +
                     sqrt(3.0/20.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(3,x,y,z,t,b);
                     
                     Ret(2,x,y,z,t,b) =
                     sqrt(3.0/20.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(0,x,y,z,t,b) +
                     sqrt(3.0/ 5.0) * sfunc::Y_2_m2(X,Y,Z) * Psi(2,x,y,z,t,b);
                     
                     Ret(3,x,y,z,t,b) =
                     sqrt(3.0/20.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(1,x,y,z,t,b) +
                     sqrt(3.0/ 5.0) * sfunc::Y_2_m2(X,Y,Z) * Psi(3,x,y,z,t,b);
                  }
            }
   else if (src_spin1_z == 0) {
      for (      int x=0; x<l_xSize; x++)
         for (   int y=0; y<l_ySize; y++)
            for (int z=0; z<l_zSize; z++) {
               int xyz[3] = {x,y,z};
               anaHAL::convert_origin(xyz, XYZ, XYZsize, 3);
               int X = XYZ[0], Y = XYZ[1], Z = XYZ[2];
               
               for (   int t=0; t<l_tSize; t++)
                  for (int b=0; b<l_bSize; b++) {
                     Ret(0,x,y,z,t,b) =
                     sqrt(3.0/10.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(0,x,y,z,t,b) +
                     sqrt(1.0/ 5.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(2,x,y,z,t,b);
                     
                     Ret(1,x,y,z,t,b) =
                     sqrt(3.0/10.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(1,x,y,z,t,b) +
                     sqrt(1.0/ 5.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(3,x,y,z,t,b);
                     
                     Ret(2,x,y,z,t,b) =
                     sqrt(1.0/ 5.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(0,x,y,z,t,b) +
                     sqrt(3.0/10.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(2,x,y,z,t,b);
                     
                     Ret(3,x,y,z,t,b) =
                     sqrt(1.0/ 5.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(1,x,y,z,t,b) +
                     sqrt(3.0/10.0) * sfunc::Y_2_m1(X,Y,Z) * Psi(3,x,y,z,t,b);
                  }
            }
      Ret *= (-1);
   }
   else if (src_spin1_z == +1)
      for (      int x=0; x<l_xSize; x++)
         for (   int y=0; y<l_ySize; y++)
            for (int z=0; z<l_zSize; z++) {
               int xyz[3] = {x,y,z};
               anaHAL::convert_origin(xyz, XYZ, XYZsize, 3);
               int X = XYZ[0], Y = XYZ[1], Z = XYZ[2];
               
               for (   int t=0; t<l_tSize; t++)
                  for (int b=0; b<l_bSize; b++) {
                     Ret(0,x,y,z,t,b) =
                     sqrt(3.0/ 5.0) * sfunc::Y_2_p2(X,Y,Z) * Psi(0,x,y,z,t,b) +
                     sqrt(3.0/20.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(2,x,y,z,t,b);
                     
                     Ret(1,x,y,z,t,b) =
                     sqrt(3.0/ 5.0) * sfunc::Y_2_p2(X,Y,Z) * Psi(1,x,y,z,t,b) +
                     sqrt(3.0/20.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(3,x,y,z,t,b);
                     
                     Ret(2,x,y,z,t,b) =
                     sqrt(3.0/20.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(0,x,y,z,t,b) +
                     sqrt(1.0/10.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(2,x,y,z,t,b);
                     
                     Ret(3,x,y,z,t,b) =
                     sqrt(3.0/20.0) * sfunc::Y_2_p1(X,Y,Z) * Psi(1,x,y,z,t,b) +
                     sqrt(1.0/10.0) * sfunc::Y_2_0 (X,Y,Z) * Psi(3,x,y,z,t,b);
                  }
            }
   else ERROR_COMMENTS("Invalid spin number.");
   
   return Ret;
}

void Potential::calc_tensor_pot_T1(ComplexField_BASE &V_C,
                                   ComplexField_BASE &V_T,
                                   const ComplexField_BASE &Rm,
                                   const ComplexField_BASE &R0,
                                   const double m1, const double m2) {
   DEBUG_LOG
   int l_xSize = R0.get_xSIZE();
   int l_ySize = R0.get_ySIZE();
   int l_zSize = R0.get_zSIZE();
   int l_tSize = R0.get_tSIZE();
   V_C.mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
   V_T.mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
   V_C = 0.0; V_T = 0.0;
   
   for (int src_spin1_z = -1; src_spin1_z <= +1; src_spin1_z++) {
      ComplexField_BASE Rm_tmp = Rm.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE R0_tmp = R0.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE KR, SR, R0_S, R0_D, KR_S, KR_D, SR_S, SR_D;
      
      KR = Potential::get_potential_T1_nume(Rm_tmp, R0_tmp, m1, m2);
      SR = Potential::get_S12_Psi(R0_tmp);
      
      R0_S = R0_tmp.rot_proj(ROT_REP_A1);
      R0_D = Potential::get_YDstar_Psi(R0_tmp-R0_S, src_spin1_z);
      
      KR_S = KR.rot_proj(ROT_REP_A1);
      KR_D = Potential::get_YDstar_Psi(KR-KR_S, src_spin1_z);
      
      SR_S = SR.rot_proj(ROT_REP_A1);
      SR_D = Potential::get_YDstar_Psi(SR-SR_S, src_spin1_z);
      
      R0_S = R0_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      KR_S = KR_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      SR_S = SR_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      
      R0_D = (R0_D.snk_spin_proj(HH_OctOct, 1,+1) +
              R0_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      KR_D = (KR_D.snk_spin_proj(HH_OctOct, 1,+1) +
              KR_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      SR_D = (SR_D.snk_spin_proj(HH_OctOct, 1,+1) +
              SR_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      
      ComplexField_BASE Det = R0_S * SR_D - SR_S * R0_D;
      
      V_C += (SR_D * KR_S - SR_S * KR_D) / Det;
      V_T += (KR_D * R0_S - KR_S * R0_D) / Det;
   }
   V_C /= 3.0;
   V_T /= 3.0;
}
void Potential::calc_tensor_pot_T1(ComplexField_BASE &V_C,
                                   ComplexField_BASE &V_T,
                                   const ComplexField_BASE &Rm,
                                   const ComplexField_BASE &R0,
                                   const ComplexField_BASE &Rp,
                                   const double m1, const double m2) {
   DEBUG_LOG
   int l_xSize = R0.get_xSIZE();
   int l_ySize = R0.get_ySIZE();
   int l_zSize = R0.get_zSIZE();
   int l_tSize = R0.get_tSIZE();
   V_C.mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
   V_T.mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
   V_C = 0.0; V_T = 0.0;
   
   for (int src_spin1_z = -1; src_spin1_z <= +1; src_spin1_z++) {
      ComplexField_BASE Rm_tmp = Rm.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE R0_tmp = R0.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE Rp_tmp = Rp.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE KR, SR, R0_S, R0_D, KR_S, KR_D, SR_S, SR_D;
      
      KR = Potential::get_potential_T1_nume(Rm_tmp, R0_tmp, Rp_tmp, m1, m2);
      SR = Potential::get_S12_Psi(R0_tmp);
      
      R0_S = R0_tmp.rot_proj(ROT_REP_A1);
      R0_D = Potential::get_YDstar_Psi(R0_tmp-R0_S, src_spin1_z);
      
      KR_S = KR.rot_proj(ROT_REP_A1);
      KR_D = Potential::get_YDstar_Psi(KR-KR_S, src_spin1_z);
      
      SR_S = SR.rot_proj(ROT_REP_A1);
      SR_D = Potential::get_YDstar_Psi(SR-SR_S, src_spin1_z);
      
      R0_S = R0_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      KR_S = KR_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      SR_S = SR_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      
      R0_D = (R0_D.snk_spin_proj(HH_OctOct, 1,+1) +
              R0_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      KR_D = (KR_D.snk_spin_proj(HH_OctOct, 1,+1) +
              KR_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      SR_D = (SR_D.snk_spin_proj(HH_OctOct, 1,+1) +
              SR_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      
      ComplexField_BASE Det = R0_S * SR_D - SR_S * R0_D;
      
      V_C += (SR_D * KR_S - SR_S * KR_D) / Det;
      V_T += (KR_D * R0_S - KR_S * R0_D) / Det;
   }
   V_C /= 3.0;
   V_T /= 3.0;
}
void Potential::calc_tensor_pot_T2(ComplexField_BASE &V_C,
                                   ComplexField_BASE &V_T,
                                   const ComplexField_BASE &Rm,
                                   const ComplexField_BASE &R0,
                                   const ComplexField_BASE &Rp,
                                   const double m1, const double m2) {
   DEBUG_LOG
   int l_xSize = R0.get_xSIZE();
   int l_ySize = R0.get_ySIZE();
   int l_zSize = R0.get_zSIZE();
   int l_tSize = R0.get_tSIZE();
   V_C.mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
   V_T.mem_alloc(1, l_xSize, l_ySize, l_zSize, l_tSize, 1);
   V_C = 0.0; V_T = 0.0;
   
   for (int src_spin1_z = -1; src_spin1_z <= +1; src_spin1_z++) {
      ComplexField_BASE Rm_tmp = Rm.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE R0_tmp = R0.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE Rp_tmp = Rp.src_spin_proj(HH_OctOct, 1, src_spin1_z);
      ComplexField_BASE KR, SR, R0_S, R0_D, KR_S, KR_D, SR_S, SR_D;
      
      KR = Potential::get_potential_T1_nume(Rm_tmp, R0_tmp, Rp_tmp, m1, m2);
      SR = Potential::get_S12_Psi(R0_tmp);
      
      R0_S = R0_tmp.rot_proj(ROT_REP_A1);
      R0_D = Potential::get_YDstar_Psi(R0_tmp-R0_S, src_spin1_z);
      
      KR_S = KR.rot_proj(ROT_REP_A1);
      KR_D = Potential::get_YDstar_Psi(KR-KR_S, src_spin1_z);
      
      SR_S = SR.rot_proj(ROT_REP_A1);
      SR_D = Potential::get_YDstar_Psi(SR-SR_S, src_spin1_z);
      
      R0_S = R0_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      KR_S = KR_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      SR_S = SR_S.snk_spin_proj(HH_OctOct, 1, src_spin1_z);
      
      R0_D = (R0_D.snk_spin_proj(HH_OctOct, 1,+1) +
              R0_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      KR_D = (KR_D.snk_spin_proj(HH_OctOct, 1,+1) +
              KR_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      SR_D = (SR_D.snk_spin_proj(HH_OctOct, 1,+1) +
              SR_D.snk_spin_proj(HH_OctOct, 1,-1)) / 2.0;
      
      ComplexField_BASE Det = R0_S * SR_D - SR_S * R0_D;
      
      V_C += (SR_D * KR_S - SR_S * KR_D) / Det;
      V_T += (KR_D * R0_S - KR_S * R0_D) / Det;
   }
   V_C /= 3.0;
   V_T /= 3.0;
}
