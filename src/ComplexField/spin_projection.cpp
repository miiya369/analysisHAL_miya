//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Implementations of the spin projection for NBS wave functions
 * @author  Takaya Miyamoto
 * @since   Mon Oct 30 18:50:56 JST 2017
 */
//--------------------------------------------------------------------------

#include <ComplexField_Base.h>

//--------------------------------------------------------------------------
namespace {
   inline void sproj_src_OctOct(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   inline void sproj_src_DecDec(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   inline void sproj_src_DecOct(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   
   inline void sproj_snk_OctOct(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   inline void sproj_snk_DecDec(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   inline void sproj_snk_DecOct(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   
   inline void sproj_src_OctMps(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
   inline void sproj_snk_OctMps(const ComplexField_BASE&, ComplexField_BASE&,
                                int spin, int spin_z);
}
//--------------------------------------------------------------------------

ComplexField_BASE ComplexField_BASE::src_spin_proj(const int idf_src_HH,
                                                   const int spin_src,
                                                   const int spin_z_src) const {
   DEBUG_LOG
   ComplexField_BASE Ret((*this).get_aSIZE(), (*this).get_xSIZE(),
                         (*this).get_ySIZE(), (*this).get_zSIZE(),
                         (*this).get_tSIZE(),                  1);
   
   if      (idf_src_HH == HH_OctOct)
      sproj_src_OctOct((*this), Ret, spin_src, spin_z_src);
   
   else if (idf_src_HH == HH_DecDec)
      sproj_src_DecDec((*this), Ret, spin_src, spin_z_src);
   
   else if (idf_src_HH == HH_DecOct)
      sproj_src_DecOct((*this), Ret, spin_src, spin_z_src);
   
   else if (idf_src_HH == HH_OctMps)
      sproj_src_OctMps((*this), Ret, spin_src, spin_z_src);
   
   else ERROR_COMMENTS("Unkown identifier for source-HH");
   
   return Ret;
}
ComplexField_BASE ComplexField_BASE::snk_spin_proj(const int idf_snk_HH,
                                                   const int spin_snk,
                                                   const int spin_z_snk) const {
   DEBUG_LOG
   if ((*this).get_bSIZE() != 1)
      ERROR_COMMENTS("Only b=1 is allowed for sink spin projection.");
   
   ComplexField_BASE Ret(                  1, (*this).get_xSIZE(),
                         (*this).get_ySIZE(), (*this).get_zSIZE(),
                         (*this).get_tSIZE(),                  1);
   
   if      (idf_snk_HH == HH_OctOct)
      sproj_snk_OctOct((*this), Ret, spin_snk, spin_z_snk);
   
   else if (idf_snk_HH == HH_DecDec)
      sproj_snk_DecDec((*this), Ret, spin_snk, spin_z_snk);
   
   else if (idf_snk_HH == HH_DecOct)
      sproj_snk_DecOct((*this), Ret, spin_snk, spin_z_snk);
   
   else if (idf_snk_HH == HH_OctMps)
      sproj_snk_OctMps((*this), Ret, spin_snk, spin_z_snk);
   
   else ERROR_COMMENTS("Unkown identifier for sink-HH");
   
   return Ret;
}
ComplexField_BASE ComplexField_BASE::spin_proj(const int idf_snk_HH,
                                               const int spin_snk,
                                               const int spin_z_snk,
                                               const int idf_src_HH,
                                               const int spin_src,
                                               const int spin_z_src) const {
   DEBUG_LOG
   
   return ((*this).
           src_spin_proj(idf_src_HH, spin_src, spin_z_src).
           snk_spin_proj(idf_snk_HH, spin_snk, spin_z_snk));
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
namespace {
   inline void sproj_src_OctOct(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_bSIZE() != 4) ERROR_COMMENTS("Invalid #.external index");
      int nsize = prj_wave.data_size();
      
#define Win(iloop, a, b) wave_in(iloop + nsize * (a + 2*b))
#define Wpr(iloop)      prj_wave(iloop)
      // spin = 0
      if      (spin == 0 && spin_z ==  0) for (int n=0; n<nsize; n++)
            Wpr(n) = (Win(n,0,1) - Win(n,1,0)) / sqrt(2);
      
      // spin = 1
      else if (spin == 1 && spin_z == +1) for (int n=0; n<nsize; n++)
            Wpr(n) =  Win(n,0,0);
      else if (spin == 1 && spin_z ==  0) for (int n=0; n<nsize; n++)
            Wpr(n) = (Win(n,0,1) + Win(n,1,0)) / sqrt(2);
      else if (spin == 1 && spin_z == -1) for (int n=0; n<nsize; n++)
            Wpr(n) =  Win(n,1,1);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   inline void sproj_snk_OctOct(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_aSIZE() != 4) ERROR_COMMENTS("Invalid #.internal index");
      int nsize = prj_wave.data_size();
      
#define Win(a, b, iloop) wave_in(a + 2*(b + 2*iloop))
#define Wpr(iloop)      prj_wave(iloop)
      // spin = 0
      if      (spin == 0 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (Win(0,1,n) - Win(1,0,n)) / sqrt(2);
      
      // spin = 1
      else if (spin == 1 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) =  Win(0,0,n);
      else if (spin == 1 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (Win(0,1,n) + Win(1,0,n)) / sqrt(2);
      else if (spin == 1 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) =  Win(1,1,n);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   //--------------------------------------------------------------------------
   inline void sproj_src_DecOct(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_bSIZE() != 12) ERROR_COMMENTS("Invalid #.external index");
      int nsize = prj_wave.data_size();
      
#define Win(iloop, a, b, m) wave_in(iloop + nsize * (a + 2*(b + 2*m)))
#define Wpr(iloop)         prj_wave(iloop)
      // spin = 1
      if      (spin == 1 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(n,1,0,0)
                   + 3.0 * COMP_IMAG * Win(n,0,1,0)
                   + 1.0 *             Win(n,1,0,1)
                   - 3.0 *             Win(n,0,1,1)
                   + 2.0 * COMP_IMAG * Win(n,0,0,2)
                   ) / (4.0 * sqrt(3));
      else if (spin == 1 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(n,0,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,1,0)
                   + 1.0 *             Win(n,0,0,1)
                   - 1.0 *             Win(n,1,1,1)
                   + 2.0 * COMP_IMAG * Win(n,1,0,2)
                   - 2.0 * COMP_IMAG * Win(n,0,1,2)
                   ) / (2.0 * sqrt(6));
      else if (spin == 1 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 3.0 * COMP_IMAG * Win(n,1,0,0)
                   - 1.0 * COMP_IMAG * Win(n,0,1,0)
                   + 3.0 *             Win(n,1,0,1)
                   - 1.0 *             Win(n,0,1,1)
                   - 2.0 * COMP_IMAG * Win(n,1,1,2)
                   ) / (4.0 * sqrt(3));
      
      // spin = 2
      else if (spin == 2 && spin_z == +2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(n,0,0,0)
                   - 1.0 *             Win(n,0,0,1)
                   ) / (2.0);
      else if (spin == 2 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(n,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,0,1,0)
                   - 1.0 *             Win(n,1,0,1)
                   - 1.0 *             Win(n,0,1,1)
                   - 2.0 * COMP_IMAG * Win(n,0,0,2)
                   ) / (4.0);
      else if (spin == 2 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(n,0,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,1,0)
                   - 1.0 *             Win(n,0,0,1)
                   - 1.0 *             Win(n,1,1,1)
                   - 2.0 * COMP_IMAG * Win(n,1,0,2)
                   - 2.0 * COMP_IMAG * Win(n,0,1,2)
                   ) / (2.0 * sqrt(6));
      else if (spin == 2 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(n,1,0,0)
                   - 1.0 * COMP_IMAG * Win(n,0,1,0)
                   - 1.0 *             Win(n,1,0,1)
                   - 1.0 *             Win(n,0,1,1)
                   - 2.0 * COMP_IMAG * Win(n,1,1,2)
                   ) / (4.0);
      else if (spin == 2 && spin_z == -2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(n,1,1,0)
                   - 1.0 *             Win(n,1,1,1)
                   ) / (2.0);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   inline void sproj_snk_DecOct(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_aSIZE() != 12) ERROR_COMMENTS("Invalid #.internal index");
      int nsize = prj_wave.data_size();
      
#define Win(a, b, m, iloop) wave_in(a + 2*(b + 2*(m + 3*iloop)))
#define Wpr(iloop)         prj_wave(iloop)
      // spin = 1
      if      (spin == 1 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(1,0,0,n)
                   + 3.0 * COMP_IMAG * Win(0,1,0,n)
                   - 1.0 *             Win(1,0,1,n)
                   + 3.0 *             Win(0,1,1,n)
                   + 2.0 * COMP_IMAG * Win(0,0,2,n)
                   ) / (4.0 * sqrt(3));
      else if (spin == 1 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(0,0,0,n)
                   + 1.0 * COMP_IMAG * Win(1,1,0,n)
                   - 1.0 *             Win(0,0,1,n)
                   + 1.0 *             Win(1,1,1,n)
                   + 2.0 * COMP_IMAG * Win(1,0,2,n)
                   - 2.0 * COMP_IMAG * Win(0,1,2,n)
                   ) / (2.0 * sqrt(6));
      else if (spin == 1 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 3.0 * COMP_IMAG * Win(1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(0,1,0,n)
                   - 3.0 *             Win(1,0,1,n)
                   + 1.0 *             Win(0,1,1,n)
                   - 2.0 * COMP_IMAG * Win(1,1,2,n)
                   ) / (4.0 * sqrt(3));
      
      // spin = 2
      else if (spin == 2 && spin_z == +2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(0,0,0,n)
                   + 1.0 *             Win(0,0,1,n)
                   ) / (2.0);
      else if (spin == 2 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(1,0,0,n)
                   + 1.0 * COMP_IMAG * Win(0,1,0,n)
                   + 1.0 *             Win(1,0,1,n)
                   + 1.0 *             Win(0,1,1,n)
                   - 2.0 * COMP_IMAG * Win(0,0,2,n)
                   ) / (4.0);
      else if (spin == 2 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(0,0,0,n)
                   + 1.0 * COMP_IMAG * Win(1,1,0,n)
                   + 1.0 *             Win(0,0,1,n)
                   + 1.0 *             Win(1,1,1,n)
                   - 2.0 * COMP_IMAG * Win(1,0,2,n)
                   - 2.0 * COMP_IMAG * Win(0,1,2,n)
                   ) / (2.0 * sqrt(6));
      else if (spin == 2 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(0,1,0,n)
                   + 1.0 *             Win(1,0,1,n)
                   + 1.0 *             Win(0,1,1,n)
                   - 2.0 * COMP_IMAG * Win(1,1,2,n)
                   ) / (4.0);
      else if (spin == 2 && spin_z == -2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(1,1,0,n)
                   + 1.0 *             Win(1,1,1,n)
                   ) / (2.0);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   //--------------------------------------------------------------------------
   inline void sproj_src_DecDec(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_bSIZE() != 36) ERROR_COMMENTS("Invalid #.external index");
      int nsize = prj_wave.data_size();
      
#define Win(iloop, a, b, m, n) wave_in(iloop + nsize * (a + 2*(b + 2*(m + 3*n))))
#define Wpr(iloop)            prj_wave(iloop)
      // spin = 0
      if      (spin == 0 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 2.0 *             Win(n,1,0,0,0)
                   + 2.0 *             Win(n,0,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,0,1,0)
                   + 1.0 * COMP_IMAG * Win(n,0,1,1,0)
                   + 1.0 *             Win(n,0,0,2,0)
                   + 1.0 *             Win(n,1,1,2,0)
                   - 1.0 * COMP_IMAG * Win(n,1,0,0,1)
                   - 1.0 * COMP_IMAG * Win(n,0,1,0,1)
                   - 2.0 *             Win(n,1,0,1,1)
                   + 2.0 *             Win(n,0,1,1,1)
                   - 1.0 * COMP_IMAG * Win(n,0,0,2,1)
                   + 1.0 * COMP_IMAG * Win(n,1,1,2,1)
                   - 1.0 *             Win(n,0,0,0,2)
                   - 1.0 *             Win(n,1,1,0,2)
                   + 1.0 * COMP_IMAG * Win(n,0,0,1,2)
                   - 1.0 * COMP_IMAG * Win(n,1,1,1,2)
                   - 2.0 *             Win(n,1,0,2,2)
                   + 2.0 *             Win(n,0,1,2,2)
                   ) / (12.0);
      
      // spin = 1
      else if (spin == 1 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 3.0 *             Win(n,0,0,0,0) + 1.0 *             Win(n,1,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,1,1,0) + 3.0 *             Win(n,1,0,2,0)
                   - 2.0 *             Win(n,0,1,2,0) + 1.0 * COMP_IMAG * Win(n,1,1,0,1)
                   + 3.0 *             Win(n,0,0,1,1) - 1.0 *             Win(n,1,1,1,1)
                   + 3.0 * COMP_IMAG * Win(n,1,0,2,1) - 2.0 * COMP_IMAG * Win(n,0,1,2,1)
                   - 2.0 *             Win(n,1,0,0,2) + 3.0 *             Win(n,0,1,0,2)
                   - 2.0 * COMP_IMAG * Win(n,1,0,1,2) + 3.0 * COMP_IMAG * Win(n,0,1,1,2)
                   + 4.0 *             Win(n,0,0,2,2)
                   ) / (6.0 * sqrt(10));
      else if (spin == 1 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 4.0 *             Win(n,1,0,0,0) + 4.0 *             Win(n,0,1,0,0)
                   - 5.0 * COMP_IMAG * Win(n,1,0,1,0) + 5.0 * COMP_IMAG * Win(n,0,1,1,0)
                   + 1.0 *             Win(n,0,0,2,0) - 1.0 *             Win(n,1,1,2,0)
                   + 5.0 * COMP_IMAG * Win(n,1,0,0,1) - 5.0 * COMP_IMAG * Win(n,0,1,0,1)
                   + 4.0 *             Win(n,1,0,1,1) + 4.0 *             Win(n,0,1,1,1)
                   - 1.0 * COMP_IMAG * Win(n,0,0,2,1) - 1.0 * COMP_IMAG * Win(n,1,1,2,1)
                   + 1.0 *             Win(n,0,0,0,2) - 1.0 *             Win(n,1,1,0,2)
                   - 1.0 * COMP_IMAG * Win(n,0,0,1,2) - 1.0 * COMP_IMAG * Win(n,1,1,1,2)
                   + 2.0 *             Win(n,1,0,2,2) + 2.0 *             Win(n,0,1,2,2)
                   ) / (12.0 * sqrt(5));
      else if (spin == 1 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 *             Win(n,0,0,0,0) + 3.0 *             Win(n,1,1,0,0)
                   - 1.0 * COMP_IMAG * Win(n,0,0,1,0) + 2.0 *             Win(n,1,0,2,0)
                   - 3.0 *             Win(n,0,1,2,0) - 1.0 * COMP_IMAG * Win(n,0,0,0,1)
                   - 1.0 *             Win(n,0,0,1,1) + 3.0 *             Win(n,1,1,1,1)
                   - 2.0 * COMP_IMAG * Win(n,1,0,2,1) + 3.0 * COMP_IMAG * Win(n,0,1,2,1)
                   - 3.0 *             Win(n,1,0,0,2) + 2.0 *             Win(n,0,1,0,2)
                   + 3.0 * COMP_IMAG * Win(n,1,0,1,2) - 2.0 * COMP_IMAG * Win(n,0,1,1,2)
                   + 4.0 *             Win(n,1,1,2,2)
                   ) / (6.0 * sqrt(10));
      
      // spin = 2
      else if (spin == 2 && spin_z == +2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 *             Win(n,1,0,0,0) - 1.0 *             Win(n,0,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,0,1,0) - 1.0 * COMP_IMAG * Win(n,0,1,1,0)
                   - 2.0 *             Win(n,0,0,2,0) + 1.0 * COMP_IMAG * Win(n,1,0,0,1)
                   - 1.0 * COMP_IMAG * Win(n,0,1,0,1) - 1.0 *             Win(n,1,0,1,1)
                   + 1.0 *             Win(n,0,1,1,1) - 2.0 * COMP_IMAG * Win(n,0,0,2,1)
                   + 2.0 *             Win(n,0,0,0,2) + 2.0 * COMP_IMAG * Win(n,0,0,1,2)
                   ) / (4.0 * sqrt(6));
      else if (spin == 2 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(n,0,0,1,0) - 1.0 *             Win(n,1,0,2,0)
                   - 1.0 * COMP_IMAG * Win(n,0,0,0,1) - 1.0 * COMP_IMAG * Win(n,1,0,2,1)
                   + 1.0 *             Win(n,0,1,0,2) + 1.0 * COMP_IMAG * Win(n,0,1,1,2)
                   ) / (2.0 * sqrt(6));
      else if (spin == 2 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,1,0,0,0) + 1.0 *             Win(n,0,1,0,0)
                   + 2.0 * COMP_IMAG * Win(n,1,0,1,0) + 2.0 * COMP_IMAG * Win(n,0,1,1,0)
                   - 1.0 *             Win(n,0,0,2,0) - 1.0 *             Win(n,1,1,2,0)
                   - 2.0 * COMP_IMAG * Win(n,1,0,0,1) - 2.0 * COMP_IMAG * Win(n,0,1,0,1)
                   - 1.0 *             Win(n,1,0,1,1) + 1.0 *             Win(n,0,1,1,1)
                   + 1.0 * COMP_IMAG * Win(n,0,0,2,1) - 1.0 * COMP_IMAG * Win(n,1,1,2,1)
                   + 1.0 *             Win(n,0,0,0,2) + 1.0 *             Win(n,1,1,0,2)
                   - 1.0 * COMP_IMAG * Win(n,0,0,1,2) + 1.0 * COMP_IMAG * Win(n,1,1,1,2)
                   + 2.0 *             Win(n,1,0,2,2) - 2.0 *             Win(n,0,1,2,2)
                   ) / (12.0);
      else if (spin == 2 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 * COMP_IMAG * Win(n,1,1,1,0) - 1.0 *             Win(n,0,1,2,0)
                   - 1.0 * COMP_IMAG * Win(n,1,1,0,1) + 1.0 * COMP_IMAG * Win(n,0,1,2,1)
                   + 1.0 *             Win(n,1,0,0,2) - 1.0 * COMP_IMAG * Win(n,1,0,1,2)
                   ) / (2.0 * sqrt(6));
      else if (spin == 2 && spin_z == -2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,1,0,0,0) + 1.0 *             Win(n,0,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,0,1,0) - 1.0 * COMP_IMAG * Win(n,0,1,1,0)
                   + 2.0 *             Win(n,1,1,2,0) + 1.0 * COMP_IMAG * Win(n,1,0,0,1)
                   - 1.0 * COMP_IMAG * Win(n,0,1,0,1) + 1.0 *             Win(n,1,0,1,1)
                   - 1.0 *             Win(n,0,1,1,1) - 2.0 * COMP_IMAG * Win(n,1,1,2,1)
                   - 2.0 *             Win(n,1,1,0,2) + 2.0 * COMP_IMAG * Win(n,1,1,1,2)
                   ) / (4.0 * sqrt(6));
      
      // spin = 3
      else if (spin == 3 && spin_z == +3) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,0,0,0,0) - 1.0 * COMP_IMAG * Win(n,0,0,1,0)
                   - 1.0 * COMP_IMAG * Win(n,0,0,0,1) + 1.0 *             Win(n,0,0,1,1)
                   ) / (4.0);
      else if (spin == 3 && spin_z == +2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,1,0,0,0) - 1.0 *             Win(n,0,1,0,0)
                   - 1.0 * COMP_IMAG * Win(n,1,0,1,0) - 1.0 * COMP_IMAG * Win(n,0,1,1,0)
                   + 2.0 *             Win(n,0,0,2,0) - 1.0 * COMP_IMAG * Win(n,1,0,0,1)
                   - 1.0 * COMP_IMAG * Win(n,0,1,0,1) + 1.0 *             Win(n,1,0,1,1)
                   + 1.0 *             Win(n,0,1,1,1) + 2.0 * COMP_IMAG * Win(n,0,0,2,1)
                   + 2.0 *             Win(n,0,0,0,2) + 2.0 * COMP_IMAG * Win(n,0,0,1,2)
                   ) / (4.0 * sqrt(6));
      else if (spin == 3 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 2.0 *             Win(n,0,0,0,0) - 1.0 *             Win(n,1,1,0,0)
                   - 1.0 * COMP_IMAG * Win(n,1,1,1,0) + 2.0 *             Win(n,1,0,2,0)
                   + 2.0 *             Win(n,0,1,2,0) - 1.0 * COMP_IMAG * Win(n,1,1,0,1)
                   + 2.0 *             Win(n,0,0,1,1) + 1.0 *             Win(n,1,1,1,1)
                   + 2.0 * COMP_IMAG * Win(n,1,0,2,1) + 2.0 * COMP_IMAG * Win(n,0,1,2,1)
                   + 2.0 *             Win(n,1,0,0,2) + 2.0 *             Win(n,0,1,0,2)
                   + 2.0 * COMP_IMAG * Win(n,1,0,1,2) + 2.0 * COMP_IMAG * Win(n,0,1,1,2)
                   - 4.0 *             Win(n,0,0,2,2)
                   ) / (4.0 * sqrt(15));
      else if (spin == 3 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 *             Win(n,1,0,0,0) + 1.0 *             Win(n,0,1,0,0)
                   - 1.0 *             Win(n,0,0,2,0) + 1.0 *             Win(n,1,1,2,0)
                   + 1.0 *             Win(n,1,0,1,1) + 1.0 *             Win(n,0,1,1,1)
                   + 1.0 * COMP_IMAG * Win(n,0,0,2,1) + 1.0 * COMP_IMAG * Win(n,1,1,2,1)
                   - 1.0 *             Win(n,0,0,0,2) + 1.0 *             Win(n,1,1,0,2)
                   + 1.0 * COMP_IMAG * Win(n,0,0,1,2) + 1.0 * COMP_IMAG * Win(n,1,1,1,2)
                   - 2.0 *             Win(n,1,0,2,2) - 2.0 *             Win(n,0,1,2,2)
                   ) / (4.0);
      else if (spin == 3 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,0,0,0,0) + 2.0 *             Win(n,1,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,0,0,1,0) - 2.0 *             Win(n,1,0,2,0)
                   - 2.0 *             Win(n,0,1,2,0) + 1.0 * COMP_IMAG * Win(n,0,0,0,1)
                   + 1.0 *             Win(n,0,0,1,1) + 2.0 *             Win(n,1,1,1,1)
                   + 2.0 * COMP_IMAG * Win(n,1,0,2,1) + 2.0 * COMP_IMAG * Win(n,0,1,2,1)
                   - 2.0 *             Win(n,1,0,0,2) - 2.0 *             Win(n,0,1,0,2)
                   + 2.0 * COMP_IMAG * Win(n,1,0,1,2) + 2.0 * COMP_IMAG * Win(n,0,1,1,2)
                   - 4.0 *             Win(n,1,1,2,2)
                   ) / (4.0 * sqrt(15));
      else if (spin == 3 && spin_z == -2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,1,0,0,0) - 1.0 *             Win(n,0,1,0,0)
                   + 1.0 * COMP_IMAG * Win(n,1,0,1,0) + 1.0 * COMP_IMAG * Win(n,0,1,1,0)
                   - 2.0 *             Win(n,1,1,2,0) + 1.0 * COMP_IMAG * Win(n,1,0,0,1)
                   + 1.0 * COMP_IMAG * Win(n,0,1,0,1) + 1.0 *             Win(n,1,0,1,1)
                   + 1.0 *             Win(n,0,1,1,1) + 2.0 * COMP_IMAG * Win(n,1,1,2,1)
                   - 2.0 *             Win(n,1,1,0,2) + 2.0 * COMP_IMAG * Win(n,1,1,1,2)
                   ) / (4.0 * sqrt(6));
      else if (spin == 3 && spin_z == -3) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(n,1,1,0,0) + 1.0 * COMP_IMAG * Win(n,1,1,1,0)
                   + 1.0 * COMP_IMAG * Win(n,1,1,0,1) + 1.0 *             Win(n,1,1,1,1)
                   ) / (4.0);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   inline void sproj_snk_DecDec(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_aSIZE() != 36) ERROR_COMMENTS("Invalid #.internal index");
      int nsize = prj_wave.data_size();
      
#define Win(a, b, m, n, iloop) wave_in(a + 2*(b + 2*(m + 3*(n + 3*iloop))))
#define Wpr(iloop)            prj_wave(iloop)
      // spin = 0
      if      (spin == 0 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 2.0 *             Win(1,0,0,0,n) + 2.0 *             Win(0,1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(1,0,1,0,n) - 1.0 * COMP_IMAG * Win(0,1,1,0,n)
                   + 1.0 *             Win(0,0,2,0,n) + 1.0 *             Win(1,1,2,0,n)
                   + 1.0 * COMP_IMAG * Win(1,0,0,1,n) + 1.0 * COMP_IMAG * Win(0,1,0,1,n)
                   - 2.0 *             Win(1,0,1,1,n) + 2.0 *             Win(0,1,1,1,n)
                   + 1.0 * COMP_IMAG * Win(0,0,2,1,n) - 1.0 * COMP_IMAG * Win(1,1,2,1,n)
                   - 1.0 *             Win(0,0,0,2,n) - 1.0 *             Win(1,1,0,2,n)
                   - 1.0 * COMP_IMAG * Win(0,0,1,2,n) + 1.0 * COMP_IMAG * Win(1,1,1,2,n)
                   - 2.0 *             Win(1,0,2,2,n) + 2.0 *             Win(0,1,2,2,n)
                   ) / (12.0);
      
      // spin = 1
      else if (spin == 1 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 3.0 *             Win(0,0,0,0,n) + 1.0 *             Win(1,1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(1,1,1,0,n) + 3.0 *             Win(1,0,2,0,n)
                   - 2.0 *             Win(0,1,2,0,n) - 1.0 * COMP_IMAG * Win(1,1,0,1,n)
                   + 3.0 *             Win(0,0,1,1,n) - 1.0 *             Win(1,1,1,1,n)
                   - 3.0 * COMP_IMAG * Win(1,0,2,1,n) + 2.0 * COMP_IMAG * Win(0,1,2,1,n)
                   - 2.0 *             Win(1,0,0,2,n) + 3.0 *             Win(0,1,0,2,n)
                   + 2.0 * COMP_IMAG * Win(1,0,1,2,n) - 3.0 * COMP_IMAG * Win(0,1,1,2,n)
                   + 4.0 *             Win(0,0,2,2,n)
                   ) / (6.0 * sqrt(10));
      else if (spin == 1 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 4.0 *             Win(1,0,0,0,n) + 4.0 *             Win(0,1,0,0,n)
                   + 5.0 * COMP_IMAG * Win(1,0,1,0,n) - 5.0 * COMP_IMAG * Win(0,1,1,0,n)
                   + 1.0 *             Win(0,0,2,0,n) - 1.0 *             Win(1,1,2,0,n)
                   - 5.0 * COMP_IMAG * Win(1,0,0,1,n) + 5.0 * COMP_IMAG * Win(0,1,0,1,n)
                   + 4.0 *             Win(1,0,1,1,n) + 4.0 *             Win(0,1,1,1,n)
                   + 1.0 * COMP_IMAG * Win(0,0,2,1,n) + 1.0 * COMP_IMAG * Win(1,1,2,1,n)
                   + 1.0 *             Win(0,0,0,2,n) - 1.0 *             Win(1,1,0,2,n)
                   + 1.0 * COMP_IMAG * Win(0,0,1,2,n) + 1.0 * COMP_IMAG * Win(1,1,1,2,n)
                   + 2.0 *             Win(1,0,2,2,n) + 2.0 *             Win(0,1,2,2,n)
                   ) / (12.0 * sqrt(5));
      else if (spin == 1 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 *             Win(0,0,0,0,n) + 3.0 *             Win(1,1,0,0,n)
                   + 1.0 * COMP_IMAG * Win(0,0,1,0,n) + 2.0 *             Win(1,0,2,0,n)
                   - 3.0 *             Win(0,1,2,0,n) + 1.0 * COMP_IMAG * Win(0,0,0,1,n)
                   - 1.0 *             Win(0,0,1,1,n) + 3.0 *             Win(1,1,1,1,n)
                   + 2.0 * COMP_IMAG * Win(1,0,2,1,n) - 3.0 * COMP_IMAG * Win(0,1,2,1,n)
                   - 3.0 *             Win(1,0,0,2,n) + 2.0 *             Win(0,1,0,2,n)
                   - 3.0 * COMP_IMAG * Win(1,0,1,2,n) + 2.0 * COMP_IMAG * Win(0,1,1,2,n)
                   + 4.0 *             Win(1,1,2,2,n)
                   ) / (6.0 * sqrt(10));
      
      // spin = 2
      else if (spin == 2 && spin_z == +2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 *             Win(1,0,0,0,n) - 1.0 *             Win(0,1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(1,0,1,0,n) + 1.0 * COMP_IMAG * Win(0,1,1,0,n)
                   - 2.0 *             Win(0,0,2,0,n) - 1.0 * COMP_IMAG * Win(1,0,0,1,n)
                   + 1.0 * COMP_IMAG * Win(0,1,0,1,n) - 1.0 *             Win(1,0,1,1,n)
                   + 1.0 *             Win(0,1,1,1,n) + 2.0 * COMP_IMAG * Win(0,0,2,1,n)
                   + 2.0 *             Win(0,0,0,2,n) - 2.0 * COMP_IMAG * Win(0,0,1,2,n)
                   ) / (4.0 * sqrt(6));
      else if (spin == 2 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(0,0,1,0,n) - 1.0 *             Win(1,0,2,0,n)
                   + 1.0 * COMP_IMAG * Win(0,0,0,1,n) + 1.0 * COMP_IMAG * Win(1,0,2,1,n)
                   + 1.0 *             Win(0,1,0,2,n) - 1.0 * COMP_IMAG * Win(0,1,1,2,n)
                   ) / (2.0 * sqrt(6));
      else if (spin == 2 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(1,0,0,0,n) + 1.0 *             Win(0,1,0,0,n)
                   - 2.0 * COMP_IMAG * Win(1,0,1,0,n) - 2.0 * COMP_IMAG * Win(0,1,1,0,n)
                   - 1.0 *             Win(0,0,2,0,n) - 1.0 *             Win(1,1,2,0,n)
                   + 2.0 * COMP_IMAG * Win(1,0,0,1,n) + 2.0 * COMP_IMAG * Win(0,1,0,1,n)
                   - 1.0 *             Win(1,0,1,1,n) + 1.0 *             Win(0,1,1,1,n)
                   - 1.0 * COMP_IMAG * Win(0,0,2,1,n) + 1.0 * COMP_IMAG * Win(1,1,2,1,n)
                   + 1.0 *             Win(0,0,0,2,n) + 1.0 *             Win(1,1,0,2,n)
                   + 1.0 * COMP_IMAG * Win(0,0,1,2,n) - 1.0 * COMP_IMAG * Win(1,1,1,2,n)
                   + 2.0 *             Win(1,0,2,2,n) - 2.0 *             Win(0,1,2,2,n)
                   ) / (12.0);
      else if (spin == 2 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 * COMP_IMAG * Win(1,1,1,0,n) - 1.0 *             Win(0,1,2,0,n)
                   + 1.0 * COMP_IMAG * Win(1,1,0,1,n) - 1.0 * COMP_IMAG * Win(0,1,2,1,n)
                   + 1.0 *             Win(1,0,0,2,n) + 1.0 * COMP_IMAG * Win(1,0,1,2,n)
                   ) / (2.0 * sqrt(6));
      else if (spin == 2 && spin_z == -2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(1,0,0,0,n) + 1.0 *             Win(0,1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(1,0,1,0,n) + 1.0 * COMP_IMAG * Win(0,1,1,0,n)
                   + 2.0 *             Win(1,1,2,0,n) - 1.0 * COMP_IMAG * Win(1,0,0,1,n)
                   + 1.0 * COMP_IMAG * Win(0,1,0,1,n) + 1.0 *             Win(1,0,1,1,n)
                   - 1.0 *             Win(0,1,1,1,n) + 2.0 * COMP_IMAG * Win(1,1,2,1,n)
                   - 2.0 *             Win(1,1,0,2,n) - 2.0 * COMP_IMAG * Win(1,1,1,2,n)
                   ) / (4.0 * sqrt(6));
      
      // spin = 3
      else if (spin == 3 && spin_z == +3) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(0,0,0,0,n) + 1.0 * COMP_IMAG * Win(0,0,1,0,n)
                   + 1.0 * COMP_IMAG * Win(0,0,0,1,n) + 1.0 *             Win(0,0,1,1,n)
                   ) / (4.0);
      else if (spin == 3 && spin_z == +2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(1,0,0,0,n) - 1.0 *             Win(0,1,0,0,n)
                   + 1.0 * COMP_IMAG * Win(1,0,1,0,n) + 1.0 * COMP_IMAG * Win(0,1,1,0,n)
                   + 2.0 *             Win(0,0,2,0,n) + 1.0 * COMP_IMAG * Win(1,0,0,1,n)
                   + 1.0 * COMP_IMAG * Win(0,1,0,1,n) + 1.0 *             Win(1,0,1,1,n)
                   + 1.0 *             Win(0,1,1,1,n) - 2.0 * COMP_IMAG * Win(0,0,2,1,n)
                   + 2.0 *             Win(0,0,0,2,n) - 2.0 * COMP_IMAG * Win(0,0,1,2,n)
                   ) / (4.0 * sqrt(6));
      else if (spin == 3 && spin_z == +1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 2.0 *             Win(0,0,0,0,n) - 1.0 *             Win(1,1,0,0,n)
                   + 1.0 * COMP_IMAG * Win(1,1,1,0,n) + 2.0 *             Win(1,0,2,0,n)
                   + 2.0 *             Win(0,1,2,0,n) + 1.0 * COMP_IMAG * Win(1,1,0,1,n)
                   + 2.0 *             Win(0,0,1,1,n) + 1.0 *             Win(1,1,1,1,n)
                   - 2.0 * COMP_IMAG * Win(1,0,2,1,n) - 2.0 * COMP_IMAG * Win(0,1,2,1,n)
                   + 2.0 *             Win(1,0,0,2,n) + 2.0 *             Win(0,1,0,2,n)
                   - 2.0 * COMP_IMAG * Win(1,0,1,2,n) - 2.0 * COMP_IMAG * Win(0,1,1,2,n)
                   - 4.0 *             Win(0,0,2,2,n)
                   ) / (4.0 * sqrt(15));
      else if (spin == 3 && spin_z ==  0) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   + 1.0 *             Win(1,0,0,0,n) + 1.0 *             Win(0,1,0,0,n)
                   - 1.0 *             Win(0,0,2,0,n) + 1.0 *             Win(1,1,2,0,n)
                   + 1.0 *             Win(1,0,1,1,n) + 1.0 *             Win(0,1,1,1,n)
                   - 1.0 * COMP_IMAG * Win(0,0,2,1,n) - 1.0 * COMP_IMAG * Win(1,1,2,1,n)
                   - 1.0 *             Win(0,0,0,2,n) + 1.0 *             Win(1,1,0,2,n)
                   - 1.0 * COMP_IMAG * Win(0,0,1,2,n) - 1.0 * COMP_IMAG * Win(1,1,1,2,n)
                   - 2.0 *             Win(1,0,2,2,n) - 2.0 *             Win(0,1,2,2,n)
                   ) / (4.0);
      else if (spin == 3 && spin_z == -1) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(0,0,0,0,n) + 2.0 *             Win(1,1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(0,0,1,0,n) - 2.0 *             Win(1,0,2,0,n)
                   - 2.0 *             Win(0,1,2,0,n) - 1.0 * COMP_IMAG * Win(0,0,0,1,n)
                   + 1.0 *             Win(0,0,1,1,n) + 2.0 *             Win(1,1,1,1,n)
                   - 2.0 * COMP_IMAG * Win(1,0,2,1,n) - 2.0 * COMP_IMAG * Win(0,1,2,1,n)
                   - 2.0 *             Win(1,0,0,2,n) - 2.0 *             Win(0,1,0,2,n)
                   - 2.0 * COMP_IMAG * Win(1,0,1,2,n) - 2.0 * COMP_IMAG * Win(0,1,1,2,n)
                   - 4.0 *             Win(1,1,2,2,n)
                   ) / (4.0 * sqrt(15));
      else if (spin == 3 && spin_z == -2) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(1,0,0,0,n) - 1.0 *             Win(0,1,0,0,n)
                   - 1.0 * COMP_IMAG * Win(1,0,1,0,n) - 1.0 * COMP_IMAG * Win(0,1,1,0,n)
                   - 2.0 *             Win(1,1,2,0,n) - 1.0 * COMP_IMAG * Win(1,0,0,1,n)
                   - 1.0 * COMP_IMAG * Win(0,1,0,1,n) + 1.0 *             Win(1,0,1,1,n)
                   + 1.0 *             Win(0,1,1,1,n) - 2.0 * COMP_IMAG * Win(1,1,2,1,n)
                   - 2.0 *             Win(1,1,0,2,n) - 2.0 * COMP_IMAG * Win(1,1,1,2,n)
                   ) / (4.0 * sqrt(6));
      else if (spin == 3 && spin_z == -3) for (int n=0; n<nsize; n++)
         Wpr(n) = (
                   - 1.0 *             Win(1,1,0,0,n) - 1.0 * COMP_IMAG * Win(1,1,1,0,n)
                   - 1.0 * COMP_IMAG * Win(1,1,0,1,n) + 1.0 *             Win(1,1,1,1,n)
                   ) / (4.0);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   //--------------------------------------------------------------------------
   inline void sproj_src_OctMps(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_bSIZE() != 2) ERROR_COMMENTS("Invalid #.external index");
      int nsize = prj_wave.data_size();
      
#define Win(iloop, a) wave_in(iloop + nsize * a)
#define Wpr(iloop)   prj_wave(iloop)
      // spin = 1/2, spin_z = + 1/2
      if      (spin == 12 && spin_z == +12)
         for (int n=0; n<nsize; n++) Wpr(n) = Win(n,0);
      
      // spin = 1/2, spin_z = - 1/2
      else if (spin == 12 && spin_z == -12)
         for (int n=0; n<nsize; n++) Wpr(n) = Win(n,1);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
   inline void sproj_snk_OctMps(const ComplexField_BASE& wave_in,
                                ComplexField_BASE& prj_wave,
                                int spin, int spin_z)
   {
      DEBUG_LOG
      if (wave_in.get_aSIZE() != 2) ERROR_COMMENTS("Invalid #.internal index");
      int nsize = prj_wave.data_size();
      
#define Win(a, iloop) wave_in(a + 2*iloop)
#define Wpr(iloop)   prj_wave(iloop)
      // spin = 1/2, spin_z = + 1/2
      if      (spin == 12 && spin_z == +12)
         for (int n=0; n<nsize; n++) Wpr(n) = Win(0,n);
      
      // spin = 1/2, spin_z = - 1/2
      else if (spin == 12 && spin_z == -12)
         for (int n=0; n<nsize; n++) Wpr(n) = Win(1,n);
      
      else ERROR_COMMENTS("Unkown spin");
#undef Win
#undef Wpr
   }
}