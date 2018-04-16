//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup NBS_base
 * @brief   implementation of the spin projection for NBS wave functions
 * @author  T. Miyamoto
 */
//--------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <complex.h>
#include <complex>
#include <omp.h>

using namespace std;

typedef double Float;
typedef complex<double> COMPLEX;
#undef  I
#define COMPLEX_I    (complex<double>(0.0,1.0))
#define COMPLEX_ZERO (complex<double>(0.0,0.0))
#define Real(x) (real(x))
#define Imag(x) (imag(x))
#define Conj(x) (conj(x))

inline int LessThanN(int n, int i)
{
   if(i < 0 || n <= i){
      printf("ERROR:[HAL_idx.h] LessThanN: i = %d >= n = %d\n",i, n);
   }
   return i;
}

#define Spinor(nd, alpha, x)  (LessThanN(nd, alpha) + nd*(x))
#define Dirac(     alpha, x)  (LessThanN( 4, alpha) +  4*(x))
#define Vector3d( alpha,  x)  (LessThanN( 3, alpha) +  3*(x))
#define Pauli(    alpha,  x)  (LessThanN( 2, alpha) +  2*(x))
#define Color(        c,  x)  (LessThanN( 3, c    ) +  3*(x))

enum {
   BB_OctOct,
   BB_DecOct,
   BB_DecDec,
   
   BB_MpsMps,
   BB_OctMps,
};

inline bool machine_is_little_endian();
inline void endian_convert(Float*, size_t);

//--------------------------------------------------------------------------

static void spin_projection (const Float* wave_in, Float* wave_out, size_t xyz_size, size_t dof_snk,
                             int idf_snk_BB, int idf_src_BB, int spin_snk, int spin_z_snk, int spin_src, int spin_z_src);

static inline void sproj_src_OctOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);
static inline void sproj_src_DecDec(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);
static inline void sproj_src_DecOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);

static inline void sproj_snk_OctOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);
static inline void sproj_snk_DecDec(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);
static inline void sproj_snk_DecOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);

static inline void sproj_src_OctMps(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);
static inline void sproj_snk_OctMps(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z);

//--------------------------------------------------------------------------

int main(const int argc, const char *argv[]) {
   if (argc != 10) {
      printf("\nusage: %s [iwave] [owave] [Lsize] [idf_snk_BB] [idf_src_BB] "
             "[spin_snk] [spin_z_snk] [spin_src] [spin_z_src]\n\n"
             "@@@ idf_[snk,src]_BB = oo, do, dd, op\n\n", argv[0]); return -1;
   }
   const char *ifname = argv[1];
   const char *ofname = argv[2];
   
   const int   Lsize  = atoi(argv[3]);
   
   const char *idf_snk_str = argv[4];
   const char *idf_src_str = argv[5];
   
   const int   spin_snk   = atoi(argv[6]);
   const int   spin_snk_z = atoi(argv[7]);
   
   const int   spin_src   = atoi(argv[8]);
   const int   spin_src_z = atoi(argv[9]);
   /*
    printf("%s %s %d %s %s %d %d %d %d\n",
    ifname, ofname, Lsize, idf_snk_str, idf_src_str,
    spin_snk, spin_snk_z, spin_src, spin_src_z);
    */
   
   int    idf_snk_BB = 999;
   size_t dof_snk;
   if        (strcmp(idf_snk_str, "oo") == 0) {
      idf_snk_BB = BB_OctOct;
      dof_snk    = 2   * 2;
   } else if (strcmp(idf_snk_str, "do") == 0) {
      idf_snk_BB = BB_DecOct;
      dof_snk    = 2   * 2*3;
   } else if (strcmp(idf_snk_str, "dd") == 0) {
      idf_snk_BB = BB_DecDec;
      dof_snk    = 2*3 * 2*3;
   } else if (strcmp(idf_snk_str, "op") == 0) {
      idf_snk_BB = BB_OctMps;
      dof_snk    = 2   * 1;
   } else { printf("Unkown idf_snk_BB = %s\n", idf_snk_str); return -1;}
   
   int    idf_src_BB = 999;
   size_t dof_src;
   if        (strcmp(idf_src_str, "oo") == 0) {
      idf_src_BB = BB_OctOct;
      dof_src    = 2   * 2;
   } else if (strcmp(idf_src_str, "do") == 0) {
      idf_src_BB = BB_DecOct;
      dof_src    = 2   * 2*3;
   } else if (strcmp(idf_src_str, "dd") == 0) {
      idf_src_BB = BB_DecDec;
      dof_src    = 2*3 * 2*3;
   } else if (strcmp(idf_src_str, "op") == 0) {
      idf_src_BB = BB_OctMps;
      dof_src    = 2   * 1;
   } else { printf("Unkown idf_src_BB = %s\n", idf_src_str); return -1;}
   
   const size_t xyz_size = Lsize*Lsize*Lsize;
   const size_t dof_wave = dof_snk * dof_src;
   
   Float *iwave = new Float [xyz_size * dof_wave * 2];
   Float *owave = new Float [xyz_size            * 2];
   printf("Read wave: size = %lu bytes\n", (sizeof(Float)*(xyz_size * dof_wave * 2)));
   ifstream ifs(ifname, ios::binary);
   if (!ifs) { printf("\nERROR while reading\n\n"); return -1;}
   ifs.read((char*)iwave, sizeof(Float)*(xyz_size * dof_wave * 2));
   ifs.close();
   if (machine_is_little_endian()) endian_convert(iwave, xyz_size * dof_wave * 2);
   
   if (spin_snk == spin_src && spin_snk_z == 999 && spin_src_z == 999) {
      printf("Spin avaraging...\n");
      Float Ave_factor = 1.0 / (2.0 * spin_snk + 1.0);
      for (size_t n=0; n<xyz_size * 2; n++) owave[n] = 0.0;
      for (int spin_z = -spin_snk; spin_z <= spin_snk; spin_z++) {
         Float *tmp_wave = new Float [xyz_size * 2];
         spin_projection(iwave, tmp_wave, xyz_size, dof_snk, idf_snk_BB, idf_src_BB, spin_snk, spin_z, spin_src, spin_z);
         for(size_t n=0; n<xyz_size * 2; n++) owave[n] += tmp_wave[n] * Ave_factor;
         delete [] tmp_wave;
      }
   } else
      spin_projection(iwave, owave, xyz_size, dof_snk, idf_snk_BB, idf_src_BB, spin_snk, spin_snk_z, spin_src, spin_src_z);
   
   if (machine_is_little_endian()) endian_convert(owave, xyz_size * 2);
   ofstream ofs(ofname, ios::binary);
   if (!ofs) { printf("\nERROR while writing\n\n"); return -1;}
   ofs.write((char*)owave, sizeof(Float)*(xyz_size * 2));
   ofs.close();
   
   return 0;
}

//--------------------------------------------------------------------------

static void spin_projection(const Float* wave_in, Float* wave_out, size_t xyz_size, size_t dof_snk,
                            int idf_snk_BB, int idf_src_BB, int spin_snk, int spin_z_snk, int spin_src, int spin_z_src)
{
   size_t   n_size   = xyz_size * dof_snk;
   
   const COMPLEX *wave_org = (COMPLEX*)wave_in;
   COMPLEX       *prj_wave = (COMPLEX*)wave_out;
   COMPLEX       *wave_tmp = new COMPLEX [n_size];
   
   // Source part projection
   if      (idf_src_BB == BB_OctOct) sproj_src_OctOct(wave_org, wave_tmp, n_size, spin_src, spin_z_src);
   else if (idf_src_BB == BB_DecDec) sproj_src_DecDec(wave_org, wave_tmp, n_size, spin_src, spin_z_src);
   else if (idf_src_BB == BB_DecOct) sproj_src_DecOct(wave_org, wave_tmp, n_size, spin_src, spin_z_src);
   else if (idf_src_BB == BB_OctMps) sproj_src_OctMps(wave_org, wave_tmp, n_size, spin_src, spin_z_src);
   else printf("Unkown idf_src_BB = %d\n", idf_src_BB);
   
   // Sink part projection
   if      (idf_snk_BB == BB_OctOct) sproj_snk_OctOct(wave_tmp, prj_wave, xyz_size, spin_snk, spin_z_snk);
   else if (idf_snk_BB == BB_DecDec) sproj_snk_DecDec(wave_tmp, prj_wave, xyz_size, spin_snk, spin_z_snk);
   else if (idf_snk_BB == BB_DecOct) sproj_snk_DecOct(wave_tmp, prj_wave, xyz_size, spin_snk, spin_z_snk);
   else if (idf_snk_BB == BB_OctMps) sproj_snk_OctMps(wave_tmp, prj_wave, xyz_size, spin_snk, spin_z_snk);
   else printf("Unkown idf_snk_BB = %d\n", idf_snk_BB);
   
   delete [] wave_tmp;
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

static inline void sproj_src_OctOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(iloop, a, b) wave_in[iloop + nsize * Pauli(a, b)]
#define Wpr(iloop)      prj_wave[iloop]
   
   // spin = 0
   if      (spin == 0 && spin_z ==  0) for (size_t n=0; n<nsize; n++) Wpr(n) = (Win(n,0,1) - Win(n,1,0)) / sqrt(2);
   // spin = 1
   else if (spin == 1 && spin_z == +1) for (size_t n=0; n<nsize; n++) Wpr(n) =  Win(n,0,0);
   else if (spin == 1 && spin_z ==  0) for (size_t n=0; n<nsize; n++) Wpr(n) = (Win(n,0,1) + Win(n,1,0)) / sqrt(2);
   else if (spin == 1 && spin_z == -1) for (size_t n=0; n<nsize; n++) Wpr(n) =  Win(n,1,1);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}
static inline void sproj_snk_OctOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(a, b, iloop) wave_in[Pauli(a, Pauli(b, iloop))]
#define Wpr(iloop)      prj_wave[iloop]
   
   // spin = 0
   if      (spin == 0 && spin_z ==  0) for (size_t n=0; n<nsize; n++) Wpr(n) = (Win(0,1,n) - Win(1,0,n)) / sqrt(2);
   // spin = 1
   else if (spin == 1 && spin_z == +1) for (size_t n=0; n<nsize; n++) Wpr(n) =  Win(0,0,n);
   else if (spin == 1 && spin_z ==  0) for (size_t n=0; n<nsize; n++) Wpr(n) = (Win(0,1,n) + Win(1,0,n)) / sqrt(2);
   else if (spin == 1 && spin_z == -1) for (size_t n=0; n<nsize; n++) Wpr(n) =  Win(1,1,n);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}

//--------------------------------------------------------------------------

static inline void sproj_src_DecOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(iloop, a, b, m) wave_in[iloop + nsize * Pauli(a, Pauli(b, m))]
#define Wpr(iloop)      prj_wave[iloop]
   
   // spin = 1
   if      (spin == 1 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(n,1,0,0) + 3.0 * COMPLEX_I * Win(n,0,1,0)
                + 1.0 *             Win(n,1,0,1) - 3.0 *             Win(n,0,1,1)
                + 2.0 * COMPLEX_I * Win(n,0,0,2)
                ) / (4.0 * sqrt(3));
   else if (spin == 1 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(n,0,0,0) + 1.0 * COMPLEX_I * Win(n,1,1,0)
                + 1.0 *             Win(n,0,0,1) - 1.0 *             Win(n,1,1,1)
                + 2.0 * COMPLEX_I * Win(n,1,0,2) - 2.0 * COMPLEX_I * Win(n,0,1,2)
                ) / (2.0 * sqrt(6));
   else if (spin == 1 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 3.0 * COMPLEX_I * Win(n,1,0,0) - 1.0 * COMPLEX_I * Win(n,0,1,0)
                + 3.0 *             Win(n,1,0,1) - 1.0 *             Win(n,0,1,1)
                - 2.0 * COMPLEX_I * Win(n,1,1,2)
                ) / (4.0 * sqrt(3));
   // spin = 2
   else if (spin == 2 && spin_z == +2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(n,0,0,0) - 1.0 *             Win(n,0,0,1)
                ) / (2.0);
   else if (spin == 2 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(n,1,0,0) + 1.0 * COMPLEX_I * Win(n,0,1,0)
                - 1.0 *             Win(n,1,0,1) - 1.0 *             Win(n,0,1,1)
                - 2.0 * COMPLEX_I * Win(n,0,0,2)
                ) / (4.0);
   else if (spin == 2 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(n,0,0,0) + 1.0 * COMPLEX_I * Win(n,1,1,0)
                - 1.0 *             Win(n,0,0,1) - 1.0 *             Win(n,1,1,1)
                - 2.0 * COMPLEX_I * Win(n,1,0,2) - 2.0 * COMPLEX_I * Win(n,0,1,2)
                ) / (2.0 * sqrt(6));
   else if (spin == 2 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(n,1,0,0) - 1.0 * COMPLEX_I * Win(n,0,1,0)
                - 1.0 *             Win(n,1,0,1) - 1.0 *             Win(n,0,1,1)
                - 2.0 * COMPLEX_I * Win(n,1,1,2)
                ) / (4.0);
   else if (spin == 2 && spin_z == -2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(n,1,1,0) - 1.0 *             Win(n,1,1,1)
                ) / (2.0);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}
static inline void sproj_snk_DecOct(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(a, b, m, iloop) wave_in[Pauli(a, Pauli(b, Vector3d(m, iloop)))]
#define Wpr(iloop)         prj_wave[iloop]
   
   // spin = 1
   if      (spin == 1 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(1,0,0,n) + 3.0 * COMPLEX_I * Win(0,1,0,n)
                - 1.0 *             Win(1,0,1,n) + 3.0 *             Win(0,1,1,n)
                + 2.0 * COMPLEX_I * Win(0,0,2,n)
                ) / (4.0 * sqrt(3));
   else if (spin == 1 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(0,0,0,n) + 1.0 * COMPLEX_I * Win(1,1,0,n)
                - 1.0 *             Win(0,0,1,n) + 1.0 *             Win(1,1,1,n)
                + 2.0 * COMPLEX_I * Win(1,0,2,n) - 2.0 * COMPLEX_I * Win(0,1,2,n)
                ) / (2.0 * sqrt(6));
   else if (spin == 1 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 3.0 * COMPLEX_I * Win(1,0,0,n) - 1.0 * COMPLEX_I * Win(0,1,0,n)
                - 3.0 *             Win(1,0,1,n) + 1.0 *             Win(0,1,1,n)
                - 2.0 * COMPLEX_I * Win(1,1,2,n)
                ) / (4.0 * sqrt(3));
   // spin = 2
   else if (spin == 2 && spin_z == +2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(0,0,0,n) + 1.0 *             Win(0,0,1,n)
                ) / (2.0);
   else if (spin == 2 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(1,0,0,n) + 1.0 * COMPLEX_I * Win(0,1,0,n)
                + 1.0 *             Win(1,0,1,n) + 1.0 *             Win(0,1,1,n)
                - 2.0 * COMPLEX_I * Win(0,0,2,n)
                ) / (4.0);
   else if (spin == 2 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(0,0,0,n) + 1.0 * COMPLEX_I * Win(1,1,0,n)
                + 1.0 *             Win(0,0,1,n) + 1.0 *             Win(1,1,1,n)
                - 2.0 * COMPLEX_I * Win(1,0,2,n) - 2.0 * COMPLEX_I * Win(0,1,2,n)
                ) / (2.0 * sqrt(6));
   else if (spin == 2 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(1,0,0,n) - 1.0 * COMPLEX_I * Win(0,1,0,n)
                + 1.0 *             Win(1,0,1,n) + 1.0 *             Win(0,1,1,n)
                - 2.0 * COMPLEX_I * Win(1,1,2,n)
                ) / (4.0);
   else if (spin == 2 && spin_z == -2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(1,1,0,n) + 1.0 *             Win(1,1,1,n)
                ) / (2.0);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}

//--------------------------------------------------------------------------

static inline void sproj_src_DecDec(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(iloop, a, b, m, n) wave_in[iloop + nsize * Pauli(a, Pauli(b, Vector3d(m, n)))]
#define Wpr(iloop)            prj_wave[iloop]
   
   // spin = 0
   if      (spin == 0 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 2.0 *             Win(n,1,0,0,0) + 2.0 *             Win(n,0,1,0,0)
                + 1.0 * COMPLEX_I * Win(n,1,0,1,0) + 1.0 * COMPLEX_I * Win(n,0,1,1,0)
                + 1.0 *             Win(n,0,0,2,0) + 1.0 *             Win(n,1,1,2,0)
                - 1.0 * COMPLEX_I * Win(n,1,0,0,1) - 1.0 * COMPLEX_I * Win(n,0,1,0,1)
                - 2.0 *             Win(n,1,0,1,1) + 2.0 *             Win(n,0,1,1,1)
                - 1.0 * COMPLEX_I * Win(n,0,0,2,1) + 1.0 * COMPLEX_I * Win(n,1,1,2,1)
                - 1.0 *             Win(n,0,0,0,2) - 1.0 *             Win(n,1,1,0,2)
                + 1.0 * COMPLEX_I * Win(n,0,0,1,2) - 1.0 * COMPLEX_I * Win(n,1,1,1,2)
                - 2.0 *             Win(n,1,0,2,2) + 2.0 *             Win(n,0,1,2,2)
                ) / (12.0);
   // spin = 1
   else if (spin == 1 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 3.0 *             Win(n,0,0,0,0) + 1.0 *             Win(n,1,1,0,0)
                + 1.0 * COMPLEX_I * Win(n,1,1,1,0) + 3.0 *             Win(n,1,0,2,0)
                - 2.0 *             Win(n,0,1,2,0) + 1.0 * COMPLEX_I * Win(n,1,1,0,1)
                + 3.0 *             Win(n,0,0,1,1) - 1.0 *             Win(n,1,1,1,1)
                + 3.0 * COMPLEX_I * Win(n,1,0,2,1) - 2.0 * COMPLEX_I * Win(n,0,1,2,1)
                - 2.0 *             Win(n,1,0,0,2) + 3.0 *             Win(n,0,1,0,2)
                - 2.0 * COMPLEX_I * Win(n,1,0,1,2) + 3.0 * COMPLEX_I * Win(n,0,1,1,2)
                + 4.0 *             Win(n,0,0,2,2)
                ) / (6.0 * sqrt(10));
   else if (spin == 1 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 4.0 *             Win(n,1,0,0,0) + 4.0 *             Win(n,0,1,0,0)
                - 5.0 * COMPLEX_I * Win(n,1,0,1,0) + 5.0 * COMPLEX_I * Win(n,0,1,1,0)
                + 1.0 *             Win(n,0,0,2,0) - 1.0 *             Win(n,1,1,2,0)
                + 5.0 * COMPLEX_I * Win(n,1,0,0,1) - 5.0 * COMPLEX_I * Win(n,0,1,0,1)
                + 4.0 *             Win(n,1,0,1,1) + 4.0 *             Win(n,0,1,1,1)
                - 1.0 * COMPLEX_I * Win(n,0,0,2,1) - 1.0 * COMPLEX_I * Win(n,1,1,2,1)
                + 1.0 *             Win(n,0,0,0,2) - 1.0 *             Win(n,1,1,0,2)
                - 1.0 * COMPLEX_I * Win(n,0,0,1,2) - 1.0 * COMPLEX_I * Win(n,1,1,1,2)
                + 2.0 *             Win(n,1,0,2,2) + 2.0 *             Win(n,0,1,2,2)
                ) / (12.0 * sqrt(5));
   else if (spin == 1 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 *             Win(n,0,0,0,0) + 3.0 *             Win(n,1,1,0,0)
                - 1.0 * COMPLEX_I * Win(n,0,0,1,0) + 2.0 *             Win(n,1,0,2,0)
                - 3.0 *             Win(n,0,1,2,0) - 1.0 * COMPLEX_I * Win(n,0,0,0,1)
                - 1.0 *             Win(n,0,0,1,1) + 3.0 *             Win(n,1,1,1,1)
                - 2.0 * COMPLEX_I * Win(n,1,0,2,1) + 3.0 * COMPLEX_I * Win(n,0,1,2,1)
                - 3.0 *             Win(n,1,0,0,2) + 2.0 *             Win(n,0,1,0,2)
                + 3.0 * COMPLEX_I * Win(n,1,0,1,2) - 2.0 * COMPLEX_I * Win(n,0,1,1,2)
                + 4.0 *             Win(n,1,1,2,2)
                ) / (6.0 * sqrt(10));
   // spin = 2
   else if (spin == 2 && spin_z == +2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 *             Win(n,1,0,0,0) - 1.0 *             Win(n,0,1,0,0)
                + 1.0 * COMPLEX_I * Win(n,1,0,1,0) - 1.0 * COMPLEX_I * Win(n,0,1,1,0)
                - 2.0 *             Win(n,0,0,2,0) + 1.0 * COMPLEX_I * Win(n,1,0,0,1)
                - 1.0 * COMPLEX_I * Win(n,0,1,0,1) - 1.0 *             Win(n,1,0,1,1)
                + 1.0 *             Win(n,0,1,1,1) - 2.0 * COMPLEX_I * Win(n,0,0,2,1)
                + 2.0 *             Win(n,0,0,0,2) + 2.0 * COMPLEX_I * Win(n,0,0,1,2)
                ) / (4.0 * sqrt(6));
   else if (spin == 2 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(n,0,0,1,0) - 1.0 *             Win(n,1,0,2,0)
                - 1.0 * COMPLEX_I * Win(n,0,0,0,1) - 1.0 * COMPLEX_I * Win(n,1,0,2,1)
                + 1.0 *             Win(n,0,1,0,2) + 1.0 * COMPLEX_I * Win(n,0,1,1,2)
                ) / (2.0 * sqrt(6));
   else if (spin == 2 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,1,0,0,0) + 1.0 *             Win(n,0,1,0,0)
                + 2.0 * COMPLEX_I * Win(n,1,0,1,0) + 2.0 * COMPLEX_I * Win(n,0,1,1,0)
                - 1.0 *             Win(n,0,0,2,0) - 1.0 *             Win(n,1,1,2,0)
                - 2.0 * COMPLEX_I * Win(n,1,0,0,1) - 2.0 * COMPLEX_I * Win(n,0,1,0,1)
                - 1.0 *             Win(n,1,0,1,1) + 1.0 *             Win(n,0,1,1,1)
                + 1.0 * COMPLEX_I * Win(n,0,0,2,1) - 1.0 * COMPLEX_I * Win(n,1,1,2,1)
                + 1.0 *             Win(n,0,0,0,2) + 1.0 *             Win(n,1,1,0,2)
                - 1.0 * COMPLEX_I * Win(n,0,0,1,2) + 1.0 * COMPLEX_I * Win(n,1,1,1,2)
                + 2.0 *             Win(n,1,0,2,2) - 2.0 *             Win(n,0,1,2,2)
                ) / (12.0);
   else if (spin == 2 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 * COMPLEX_I * Win(n,1,1,1,0) - 1.0 *             Win(n,0,1,2,0)
                - 1.0 * COMPLEX_I * Win(n,1,1,0,1) + 1.0 * COMPLEX_I * Win(n,0,1,2,1)
                + 1.0 *             Win(n,1,0,0,2) - 1.0 * COMPLEX_I * Win(n,1,0,1,2)
                ) / (2.0 * sqrt(6));
   else if (spin == 2 && spin_z == -2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,1,0,0,0) + 1.0 *             Win(n,0,1,0,0)
                + 1.0 * COMPLEX_I * Win(n,1,0,1,0) - 1.0 * COMPLEX_I * Win(n,0,1,1,0)
                + 2.0 *             Win(n,1,1,2,0) + 1.0 * COMPLEX_I * Win(n,1,0,0,1)
                - 1.0 * COMPLEX_I * Win(n,0,1,0,1) + 1.0 *             Win(n,1,0,1,1)
                - 1.0 *             Win(n,0,1,1,1) - 2.0 * COMPLEX_I * Win(n,1,1,2,1)
                - 2.0 *             Win(n,1,1,0,2) + 2.0 * COMPLEX_I * Win(n,1,1,1,2)
                ) / (4.0 * sqrt(6));
   // spin = 3
   else if (spin == 3 && spin_z == +3) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,0,0,0,0) - 1.0 * COMPLEX_I * Win(n,0,0,1,0)
                - 1.0 * COMPLEX_I * Win(n,0,0,0,1) + 1.0 *             Win(n,0,0,1,1)
                ) / (4.0);
   else if (spin == 3 && spin_z == +2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,1,0,0,0) - 1.0 *             Win(n,0,1,0,0)
                - 1.0 * COMPLEX_I * Win(n,1,0,1,0) - 1.0 * COMPLEX_I * Win(n,0,1,1,0)
                + 2.0 *             Win(n,0,0,2,0) - 1.0 * COMPLEX_I * Win(n,1,0,0,1)
                - 1.0 * COMPLEX_I * Win(n,0,1,0,1) + 1.0 *             Win(n,1,0,1,1)
                + 1.0 *             Win(n,0,1,1,1) + 2.0 * COMPLEX_I * Win(n,0,0,2,1)
                + 2.0 *             Win(n,0,0,0,2) + 2.0 * COMPLEX_I * Win(n,0,0,1,2)
                ) / (4.0 * sqrt(6));
   else if (spin == 3 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 2.0 *             Win(n,0,0,0,0) - 1.0 *             Win(n,1,1,0,0)
                - 1.0 * COMPLEX_I * Win(n,1,1,1,0) + 2.0 *             Win(n,1,0,2,0)
                + 2.0 *             Win(n,0,1,2,0) - 1.0 * COMPLEX_I * Win(n,1,1,0,1)
                + 2.0 *             Win(n,0,0,1,1) + 1.0 *             Win(n,1,1,1,1)
                + 2.0 * COMPLEX_I * Win(n,1,0,2,1) + 2.0 * COMPLEX_I * Win(n,0,1,2,1)
                + 2.0 *             Win(n,1,0,0,2) + 2.0 *             Win(n,0,1,0,2)
                + 2.0 * COMPLEX_I * Win(n,1,0,1,2) + 2.0 * COMPLEX_I * Win(n,0,1,1,2)
                - 4.0 *             Win(n,0,0,2,2)
                ) / (4.0 * sqrt(15));
   else if (spin == 3 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 *             Win(n,1,0,0,0) + 1.0 *             Win(n,0,1,0,0)
                - 1.0 *             Win(n,0,0,2,0) + 1.0 *             Win(n,1,1,2,0)
                + 1.0 *             Win(n,1,0,1,1) + 1.0 *             Win(n,0,1,1,1)
                + 1.0 * COMPLEX_I * Win(n,0,0,2,1) + 1.0 * COMPLEX_I * Win(n,1,1,2,1)
                - 1.0 *             Win(n,0,0,0,2) + 1.0 *             Win(n,1,1,0,2)
                + 1.0 * COMPLEX_I * Win(n,0,0,1,2) + 1.0 * COMPLEX_I * Win(n,1,1,1,2)
                - 2.0 *             Win(n,1,0,2,2) - 2.0 *             Win(n,0,1,2,2)
                ) / (4.0);
   else if (spin == 3 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,0,0,0,0) + 2.0 *             Win(n,1,1,0,0)
                + 1.0 * COMPLEX_I * Win(n,0,0,1,0) - 2.0 *             Win(n,1,0,2,0)
                - 2.0 *             Win(n,0,1,2,0) + 1.0 * COMPLEX_I * Win(n,0,0,0,1)
                + 1.0 *             Win(n,0,0,1,1) + 2.0 *             Win(n,1,1,1,1)
                + 2.0 * COMPLEX_I * Win(n,1,0,2,1) + 2.0 * COMPLEX_I * Win(n,0,1,2,1)
                - 2.0 *             Win(n,1,0,0,2) - 2.0 *             Win(n,0,1,0,2)
                + 2.0 * COMPLEX_I * Win(n,1,0,1,2) + 2.0 * COMPLEX_I * Win(n,0,1,1,2)
                - 4.0 *             Win(n,1,1,2,2)
                ) / (4.0 * sqrt(15));
   else if (spin == 3 && spin_z == -2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,1,0,0,0) - 1.0 *             Win(n,0,1,0,0)
                + 1.0 * COMPLEX_I * Win(n,1,0,1,0) + 1.0 * COMPLEX_I * Win(n,0,1,1,0)
                - 2.0 *             Win(n,1,1,2,0) + 1.0 * COMPLEX_I * Win(n,1,0,0,1)
                + 1.0 * COMPLEX_I * Win(n,0,1,0,1) + 1.0 *             Win(n,1,0,1,1)
                + 1.0 *             Win(n,0,1,1,1) + 2.0 * COMPLEX_I * Win(n,1,1,2,1)
                - 2.0 *             Win(n,1,1,0,2) + 2.0 * COMPLEX_I * Win(n,1,1,1,2)
                ) / (4.0 * sqrt(6));
   else if (spin == 3 && spin_z == -3) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(n,1,1,0,0) + 1.0 * COMPLEX_I * Win(n,1,1,1,0)
                + 1.0 * COMPLEX_I * Win(n,1,1,0,1) + 1.0 *             Win(n,1,1,1,1)
                ) / (4.0);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}
static inline void sproj_snk_DecDec(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(a, b, m, n, iloop) wave_in[Pauli(a, Pauli(b, Vector3d(m, Vector3d(n, iloop))))]
#define Wpr(iloop)            prj_wave[iloop]
   
   // spin = 0
   if      (spin == 0 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 2.0 *             Win(1,0,0,0,n) + 2.0 *             Win(0,1,0,0,n)
                - 1.0 * COMPLEX_I * Win(1,0,1,0,n) - 1.0 * COMPLEX_I * Win(0,1,1,0,n)
                + 1.0 *             Win(0,0,2,0,n) + 1.0 *             Win(1,1,2,0,n)
                + 1.0 * COMPLEX_I * Win(1,0,0,1,n) + 1.0 * COMPLEX_I * Win(0,1,0,1,n)
                - 2.0 *             Win(1,0,1,1,n) + 2.0 *             Win(0,1,1,1,n)
                + 1.0 * COMPLEX_I * Win(0,0,2,1,n) - 1.0 * COMPLEX_I * Win(1,1,2,1,n)
                - 1.0 *             Win(0,0,0,2,n) - 1.0 *             Win(1,1,0,2,n)
                - 1.0 * COMPLEX_I * Win(0,0,1,2,n) + 1.0 * COMPLEX_I * Win(1,1,1,2,n)
                - 2.0 *             Win(1,0,2,2,n) + 2.0 *             Win(0,1,2,2,n)
                ) / (12.0);
   // spin = 1
   else if (spin == 1 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 3.0 *             Win(0,0,0,0,n) + 1.0 *             Win(1,1,0,0,n)
                - 1.0 * COMPLEX_I * Win(1,1,1,0,n) + 3.0 *             Win(1,0,2,0,n)
                - 2.0 *             Win(0,1,2,0,n) - 1.0 * COMPLEX_I * Win(1,1,0,1,n)
                + 3.0 *             Win(0,0,1,1,n) - 1.0 *             Win(1,1,1,1,n)
                - 3.0 * COMPLEX_I * Win(1,0,2,1,n) + 2.0 * COMPLEX_I * Win(0,1,2,1,n)
                - 2.0 *             Win(1,0,0,2,n) + 3.0 *             Win(0,1,0,2,n)
                + 2.0 * COMPLEX_I * Win(1,0,1,2,n) - 3.0 * COMPLEX_I * Win(0,1,1,2,n)
                + 4.0 *             Win(0,0,2,2,n)
                ) / (6.0 * sqrt(10));
   else if (spin == 1 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 4.0 *             Win(1,0,0,0,n) + 4.0 *             Win(0,1,0,0,n)
                + 5.0 * COMPLEX_I * Win(1,0,1,0,n) - 5.0 * COMPLEX_I * Win(0,1,1,0,n)
                + 1.0 *             Win(0,0,2,0,n) - 1.0 *             Win(1,1,2,0,n)
                - 5.0 * COMPLEX_I * Win(1,0,0,1,n) + 5.0 * COMPLEX_I * Win(0,1,0,1,n)
                + 4.0 *             Win(1,0,1,1,n) + 4.0 *             Win(0,1,1,1,n)
                + 1.0 * COMPLEX_I * Win(0,0,2,1,n) + 1.0 * COMPLEX_I * Win(1,1,2,1,n)
                + 1.0 *             Win(0,0,0,2,n) - 1.0 *             Win(1,1,0,2,n)
                + 1.0 * COMPLEX_I * Win(0,0,1,2,n) + 1.0 * COMPLEX_I * Win(1,1,1,2,n)
                + 2.0 *             Win(1,0,2,2,n) + 2.0 *             Win(0,1,2,2,n)
                ) / (12.0 * sqrt(5));
   else if (spin == 1 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 *             Win(0,0,0,0,n) + 3.0 *             Win(1,1,0,0,n)
                + 1.0 * COMPLEX_I * Win(0,0,1,0,n) + 2.0 *             Win(1,0,2,0,n)
                - 3.0 *             Win(0,1,2,0,n) + 1.0 * COMPLEX_I * Win(0,0,0,1,n)
                - 1.0 *             Win(0,0,1,1,n) + 3.0 *             Win(1,1,1,1,n)
                + 2.0 * COMPLEX_I * Win(1,0,2,1,n) - 3.0 * COMPLEX_I * Win(0,1,2,1,n)
                - 3.0 *             Win(1,0,0,2,n) + 2.0 *             Win(0,1,0,2,n)
                - 3.0 * COMPLEX_I * Win(1,0,1,2,n) + 2.0 * COMPLEX_I * Win(0,1,1,2,n)
                + 4.0 *             Win(1,1,2,2,n)
                ) / (6.0 * sqrt(10));
   // spin = 2
   else if (spin == 2 && spin_z == +2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 *             Win(1,0,0,0,n) - 1.0 *             Win(0,1,0,0,n)
                - 1.0 * COMPLEX_I * Win(1,0,1,0,n) + 1.0 * COMPLEX_I * Win(0,1,1,0,n)
                - 2.0 *             Win(0,0,2,0,n) - 1.0 * COMPLEX_I * Win(1,0,0,1,n)
                + 1.0 * COMPLEX_I * Win(0,1,0,1,n) - 1.0 *             Win(1,0,1,1,n)
                + 1.0 *             Win(0,1,1,1,n) + 2.0 * COMPLEX_I * Win(0,0,2,1,n)
                + 2.0 *             Win(0,0,0,2,n) - 2.0 * COMPLEX_I * Win(0,0,1,2,n)
                ) / (4.0 * sqrt(6));
   else if (spin == 2 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(0,0,1,0,n) - 1.0 *             Win(1,0,2,0,n)
                + 1.0 * COMPLEX_I * Win(0,0,0,1,n) + 1.0 * COMPLEX_I * Win(1,0,2,1,n)
                + 1.0 *             Win(0,1,0,2,n) - 1.0 * COMPLEX_I * Win(0,1,1,2,n)
                ) / (2.0 * sqrt(6));
   else if (spin == 2 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(1,0,0,0,n) + 1.0 *             Win(0,1,0,0,n)
                - 2.0 * COMPLEX_I * Win(1,0,1,0,n) - 2.0 * COMPLEX_I * Win(0,1,1,0,n)
                - 1.0 *             Win(0,0,2,0,n) - 1.0 *             Win(1,1,2,0,n)
                + 2.0 * COMPLEX_I * Win(1,0,0,1,n) + 2.0 * COMPLEX_I * Win(0,1,0,1,n)
                - 1.0 *             Win(1,0,1,1,n) + 1.0 *             Win(0,1,1,1,n)
                - 1.0 * COMPLEX_I * Win(0,0,2,1,n) + 1.0 * COMPLEX_I * Win(1,1,2,1,n)
                + 1.0 *             Win(0,0,0,2,n) + 1.0 *             Win(1,1,0,2,n)
                + 1.0 * COMPLEX_I * Win(0,0,1,2,n) - 1.0 * COMPLEX_I * Win(1,1,1,2,n)
                + 2.0 *             Win(1,0,2,2,n) - 2.0 *             Win(0,1,2,2,n)
                ) / (12.0);
   else if (spin == 2 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 * COMPLEX_I * Win(1,1,1,0,n) - 1.0 *             Win(0,1,2,0,n)
                + 1.0 * COMPLEX_I * Win(1,1,0,1,n) - 1.0 * COMPLEX_I * Win(0,1,2,1,n)
                + 1.0 *             Win(1,0,0,2,n) + 1.0 * COMPLEX_I * Win(1,0,1,2,n)
                ) / (2.0 * sqrt(6));
   else if (spin == 2 && spin_z == -2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(1,0,0,0,n) + 1.0 *             Win(0,1,0,0,n)
                - 1.0 * COMPLEX_I * Win(1,0,1,0,n) + 1.0 * COMPLEX_I * Win(0,1,1,0,n)
                + 2.0 *             Win(1,1,2,0,n) - 1.0 * COMPLEX_I * Win(1,0,0,1,n)
                + 1.0 * COMPLEX_I * Win(0,1,0,1,n) + 1.0 *             Win(1,0,1,1,n)
                - 1.0 *             Win(0,1,1,1,n) + 2.0 * COMPLEX_I * Win(1,1,2,1,n)
                - 2.0 *             Win(1,1,0,2,n) - 2.0 * COMPLEX_I * Win(1,1,1,2,n)
                ) / (4.0 * sqrt(6));
   // spin = 3
   else if (spin == 3 && spin_z == +3) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(0,0,0,0,n) + 1.0 * COMPLEX_I * Win(0,0,1,0,n)
                + 1.0 * COMPLEX_I * Win(0,0,0,1,n) + 1.0 *             Win(0,0,1,1,n)
                ) / (4.0);
   else if (spin == 3 && spin_z == +2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(1,0,0,0,n) - 1.0 *             Win(0,1,0,0,n)
                + 1.0 * COMPLEX_I * Win(1,0,1,0,n) + 1.0 * COMPLEX_I * Win(0,1,1,0,n)
                + 2.0 *             Win(0,0,2,0,n) + 1.0 * COMPLEX_I * Win(1,0,0,1,n)
                + 1.0 * COMPLEX_I * Win(0,1,0,1,n) + 1.0 *             Win(1,0,1,1,n)
                + 1.0 *             Win(0,1,1,1,n) - 2.0 * COMPLEX_I * Win(0,0,2,1,n)
                + 2.0 *             Win(0,0,0,2,n) - 2.0 * COMPLEX_I * Win(0,0,1,2,n)
                ) / (4.0 * sqrt(6));
   else if (spin == 3 && spin_z == +1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 2.0 *             Win(0,0,0,0,n) - 1.0 *             Win(1,1,0,0,n)
                + 1.0 * COMPLEX_I * Win(1,1,1,0,n) + 2.0 *             Win(1,0,2,0,n)
                + 2.0 *             Win(0,1,2,0,n) + 1.0 * COMPLEX_I * Win(1,1,0,1,n)
                + 2.0 *             Win(0,0,1,1,n) + 1.0 *             Win(1,1,1,1,n)
                - 2.0 * COMPLEX_I * Win(1,0,2,1,n) - 2.0 * COMPLEX_I * Win(0,1,2,1,n)
                + 2.0 *             Win(1,0,0,2,n) + 2.0 *             Win(0,1,0,2,n)
                - 2.0 * COMPLEX_I * Win(1,0,1,2,n) - 2.0 * COMPLEX_I * Win(0,1,1,2,n)
                - 4.0 *             Win(0,0,2,2,n)
                ) / (4.0 * sqrt(15));
   else if (spin == 3 && spin_z ==  0) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                + 1.0 *             Win(1,0,0,0,n) + 1.0 *             Win(0,1,0,0,n)
                - 1.0 *             Win(0,0,2,0,n) + 1.0 *             Win(1,1,2,0,n)
                + 1.0 *             Win(1,0,1,1,n) + 1.0 *             Win(0,1,1,1,n)
                - 1.0 * COMPLEX_I * Win(0,0,2,1,n) - 1.0 * COMPLEX_I * Win(1,1,2,1,n)
                - 1.0 *             Win(0,0,0,2,n) + 1.0 *             Win(1,1,0,2,n)
                - 1.0 * COMPLEX_I * Win(0,0,1,2,n) - 1.0 * COMPLEX_I * Win(1,1,1,2,n)
                - 2.0 *             Win(1,0,2,2,n) - 2.0 *             Win(0,1,2,2,n)
                ) / (4.0);
   else if (spin == 3 && spin_z == -1) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(0,0,0,0,n) + 2.0 *             Win(1,1,0,0,n)
                - 1.0 * COMPLEX_I * Win(0,0,1,0,n) - 2.0 *             Win(1,0,2,0,n)
                - 2.0 *             Win(0,1,2,0,n) - 1.0 * COMPLEX_I * Win(0,0,0,1,n)
                + 1.0 *             Win(0,0,1,1,n) + 2.0 *             Win(1,1,1,1,n)
                - 2.0 * COMPLEX_I * Win(1,0,2,1,n) - 2.0 * COMPLEX_I * Win(0,1,2,1,n)
                - 2.0 *             Win(1,0,0,2,n) - 2.0 *             Win(0,1,0,2,n)
                - 2.0 * COMPLEX_I * Win(1,0,1,2,n) - 2.0 * COMPLEX_I * Win(0,1,1,2,n)
                - 4.0 *             Win(1,1,2,2,n)
                ) / (4.0 * sqrt(15));
   else if (spin == 3 && spin_z == -2) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(1,0,0,0,n) - 1.0 *             Win(0,1,0,0,n)
                - 1.0 * COMPLEX_I * Win(1,0,1,0,n) - 1.0 * COMPLEX_I * Win(0,1,1,0,n)
                - 2.0 *             Win(1,1,2,0,n) - 1.0 * COMPLEX_I * Win(1,0,0,1,n)
                - 1.0 * COMPLEX_I * Win(0,1,0,1,n) + 1.0 *             Win(1,0,1,1,n)
                + 1.0 *             Win(0,1,1,1,n) - 2.0 * COMPLEX_I * Win(1,1,2,1,n)
                - 2.0 *             Win(1,1,0,2,n) - 2.0 * COMPLEX_I * Win(1,1,1,2,n)
                ) / (4.0 * sqrt(6));
   else if (spin == 3 && spin_z == -3) for (size_t n=0; n<nsize; n++)
      Wpr(n) = (
                - 1.0 *             Win(1,1,0,0,n) - 1.0 * COMPLEX_I * Win(1,1,1,0,n)
                - 1.0 * COMPLEX_I * Win(1,1,0,1,n) + 1.0 *             Win(1,1,1,1,n)
                ) / (4.0);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}

//--------------------------------------------------------------------------

static inline void sproj_src_OctMps(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(iloop, a) wave_in[iloop + nsize * Pauli(a, 0)]
#define Wpr(iloop)   prj_wave[iloop]
   
   // spin = 1/2, spin_z = + 1/2
   if      (spin == 12 && spin_z == +12) for (size_t n=0; n<nsize; n++) Wpr(n) = Win(n,0);
   // spin = 1/2, spin_z = - 1/2
   else if (spin == 12 && spin_z == -12) for (size_t n=0; n<nsize; n++) Wpr(n) = Win(n,1);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}
static inline void sproj_snk_OctMps(const COMPLEX* wave_in, COMPLEX* prj_wave, size_t nsize, int spin, int spin_z)
{
#define Win(a, iloop) wave_in[Pauli(a, iloop)]
#define Wpr(iloop)   prj_wave[iloop]
   
   // spin = 1/2, spin_z = + 1/2
   if      (spin == 12 && spin_z == +12) for (size_t n=0; n<nsize; n++) Wpr(n) = Win(0,n);
   // spin = 1/2, spin_z = - 1/2
   else if (spin == 12 && spin_z == -12) for (size_t n=0; n<nsize; n++) Wpr(n) = Win(1,n);
   
   else printf("Unkown spin = %d (spin_z = %d)\n", spin, spin_z);
   
#undef Win
#undef Wpr
}

//--------------------------------------------------------------------------

bool machine_is_little_endian()
{
   int endianTEST = 1;
   if (*(char*)&endianTEST) return true;
   else                     return false;
}

void endian_convert(Float *DATA, size_t DATA_size)
{
#pragma omp parallel for
   for (size_t k=0; k<DATA_size; k++) {
      char dummy[8];
      for (int j=0; j<8; j++) dummy[j] = ((char*)&DATA[k])[j];
      for (int j=0; j<8; j++) ((char*)&DATA[k])[j] = dummy[7-j];
   }
}

