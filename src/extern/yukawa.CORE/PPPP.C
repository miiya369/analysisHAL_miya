//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup CORE_util_group
 * @brief
 * @author  N.Ishii
 * @since   Thu Aug 12 23:27:17 JST 2010
 *
 * Order of the spnor indices are <br>
 * index = alphaF + 2*(betaF + 2*(alphaI + 2*betaI))
 *
 */
//--------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/GPM.h>
#include <yukawa/matrix.h>
#include <yukawa/cubic_group.h>
#include <yukawa/misc.h>

#include <yukawa.CORE/PPPP.h>

Start_Name_Space_Yukawa

typedef std::complex<double>    Complex;
typedef yukawa::matrix<Complex> Matrix;
typedef yukawa::GPM<double>     GPMR;

static bool is_initialized = false;

static GPMR VxVxVxV;

GPMR PPPP_util::C;
GPMR PPPP_util::CxC;
GPMR PPPP_util::CxCxCxC;

Matrix PPPP_util::S[25];
Matrix PPPP_util::SxS[25];
Matrix PPPP_util::SstarxSstar[25];
Matrix PPPP_util::SxSxSstarxSstar[25];

//--------------------------------------------------------------------------
/**
 *
 */
//--------------------------------------------------------------------------

void PPPP_util::initialize()
{
  if(is_initialized == true) return;
  cubic_group::initialize();

  C       = GPMR(1,0, -1.0,1.0);
  CxC     = Kronecker_product(C, C);
  CxCxCxC = Kronecker_product(CxC, CxC);

  for(int i = 0; i < 24; i++) S[i] = Matrix(2);
  const double sqrt2 = 1.41421356237309504880;

  S[0](0,0) = 1.0;  S[0](0,1) = 0.0;
  S[0](1,0) = 0.0;  S[0](1,1) = 1.0;
  
  S[2](0,0) =  1.0/sqrt2;  S[2](0,1) = -1.0/sqrt2;
  S[2](1,0) =  1.0/sqrt2;  S[2](1,1) =  1.0/sqrt2;
  
  S[3](0,0) = Complex(1.0,-1.0)/sqrt2; S[3](0,1) = 0.0;
  S[3](1,0) = 0.0;                     S[3](1,1) = Complex(1.0,1.0)/sqrt2;
  
  cubic_group::construct_representation_matrices(S);
  
  Matrix Sstar[24];
  for(int i = 0; i < 24; i++) Sstar[i] = Matrix(2);
  for(int i = 0; i < 24; i++){
    Sstar[i](0,0) = conj(S[i](0,0));  Sstar[i](0,1) = conj(S[i](0,1));
    Sstar[i](1,0) = conj(S[i](1,0));  Sstar[i](1,1) = conj(S[i](1,1));
  }
  
  for(int i = 0; i < 24; i++){
    SxS[i]         = Kronecker_product(S[i], S[i]);
    SstarxSstar[i] = Kronecker_product(Sstar[i], Sstar[i]);
  }
  
  for(int i = 0; i < 24; i++)
    SxSxSstarxSstar[i] = Kronecker_product(SxS[i], SstarxSstar[i]);  

  S[24]               = S[0];
  SxS[24]             = SxS[0];
  SstarxSstar[24]     = SstarxSstar[0];
  SxSxSstarxSstar[24] = SxSxSstarxSstar[0];

  GPMR V(1,0, 1.0,-1.0);
  GPMR VxV = Kronecker_product(V, V);
  VxVxVxV  = Kronecker_product(VxV, VxV);

  is_initialized = true;
}

//--------------------------------------------------------------------------
/**
 * @brief spinor rotation
 *
 *
 * @f[
 *   g:
 *   G_{\alpha\beta,\alpha'\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 *   \mapsto
 *   G'_{\alpha\beta,\alpha'\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 *   \equiv
 *   S_{\alpha \tilde\alpha}(g)
 *   S_{\beta  \tilde\beta }(g)
 *   \times
 *   G_{\tilde\alpha \tilde\beta, \tilde\alpha' \tilde\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 *   \times
 *   S_{\tilde\alpha'\alpha'}(g^{-1})
 *   S_{\tilde\beta' \beta' }(g^{-1})
 * @f]
 *
 */
//--------------------------------------------------------------------------

void PPPP_util::rotate(Complex* buf, int ith)
{
  Matrix& mat(SxSxSstarxSstar[ith]);

  Complex tmp[16];
  for(int i = 0; i < 2 * 16; i++)
    ((double*)tmp)[i] = ((double*)buf)[i];
  
  for(int i = 0; i < 16; i++){
    Complex sum = 0.0;
    for(int j = 0; j < 16; j++) sum += mat(i,j) * tmp[j];
    buf[i] = sum;
  }
}

//--------------------------------------------------------------------------
/**
 * @brief charge conjugation
 *
 * multiply the charge conjugation matrices on the nucleon four point
 * correlators represented in two-component Pauli spinors as
 *
 * @f[
 *   G_{\alpha\beta;\alpha'\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 *   \mapsto
 *   C_{\alpha\tilde\alpha}
 *   C_{\beta \tilde\beta}
 *   \times
 *   G^*_{\tilde\alpha \tilde\beta; \tilde\alpha'\tilde\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 *   \times
 *   C_{\tilde\alpha'\alpha'}
 *   C_{\tilde\beta' \beta'},
 * @f]
 * with
 * @f[
 *   C
 *   \equiv
 *   \left[
 *   \begin{array}{cc}
 *   0 & -1 \\
 *   1 & 0
 *   \end{array}
 *   \right].
 * @f]
 */
//--------------------------------------------------------------------------

void PPPP_util::charge_conjugation(Complex* buf)
{
  Complex tmp[2*2 * 2*2];
  for(int i = 0; i < 16; i++)
    tmp[i] = CxCxCxC.Z[i] * conj(buf[CxCxCxC.I[i]]);

  for(int i = 0; i < 2 * 16; i++)
    ((double*)buf)[i] = ((double*)tmp)[i];
}

//---------------------------------------------------------------------
/**
 * @brief change the convention of gamma matrix from cps
 *        to the standard one. (This is releted to N.Ishii
 *        H.Nemura, and K.Murano.)
 *
 * @f$\gamma_x@f$ and @f$\gamma_z@f$ of @a CPS++ have opposite sign
 * compared to the standard ones. To convert @a CPS++ convension to
 * the standard convention, it is easy to use
 * @f[
 *   U
 *   \equiv
 *   \gamma_5\gamma_t\gamma_y
 *   =
 *   \left[
 *   \begin{array}{cc|cc}
 *    0 & 1 & 0 & 0 \\
 *   -1 & 0 & 0 & 0 \\
 *   \hline
 *    0 & 0 & 0 & 1 \\
 *    0 & 0 & -1& 0
 *   \end{array}
 *   \right].
 * @f]
 * as
 * @f[
 *   q_{\alpha}(x) \mapsto U_{\alpha\alpha'} q_{\alpha'}(x),
 *   \hspace{2em}
 *   \mbox{and}
 *   \hspace{2em}
 *   \bar q_{\alpha}(x) \mapsto U_{\alpha\alpha'} \bar q_{\alpha'}(x)
 * @f].
 * For four-point nucleon correlators, the conversion is performed as
 * @f[
 *   G_{\alpha\beta;\alpha'\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 *   \mapsto
 *   U_{\alpha \tilde\alpha}
 *   U_{\beta  \tilde\beta}
 *   U_{\alpha'\tilde\alpha'}
 *   U_{\beta' \tilde\beta'}
 *   \times
 *   G_{\tilde\alpha \tilde\beta; \tilde\alpha'\tilde\beta'}
 *   (\vec x,\vec y, t; \vec x',\vec y')
 * @f]
 */
//---------------------------------------------------------------------

void PPPP_util::cps2standard(Complex* buf)
{
  Complex tmp[16];

  for(int i = 0; i < 16; i++)
    tmp[i] = VxVxVxV.Z[i] * buf[VxVxVxV.I[i]];

    for(int i = 0; i < 2 * 16; i++)
      ((double*)buf)[i] = ((double*)tmp)[i];
}

/// read
void PPPP_util::read(Complex* a, FILE* fp, bool is_change_endian)
{
  int length = fread(a, sizeof(Complex), 2*2 * 2*2, fp);
  if(length != 2*2 * 2*2){
    fprintf(stderr,
	    "ERROR(PPPP_util::read(Complex*,FILE*,bool))\n");
    exit(1);
  }
  if(is_change_endian == true){
    yukawa::change_endian((double*)a, 2 * 2*2 * 2*2);
  }
}

/// write
void PPPP_util::write(const Complex* a, FILE* fp, bool is_change_endian)
{
  if(is_change_endian == true){
    Complex b[2*2 * 2*2];
    for(int i = 0; i < 2 * 2*2 * 2*2; i++)
      ((double*)b)[i] = ((double*)a)[i];
    yukawa::change_endian((double*)b, 2 * 2*2 * 2*2);

    int length = fwrite(b, sizeof(Complex), 2*2 * 2*2, fp);
    if(length != 2*2 * 2*2){
      fprintf(stderr,
	      "ERROR(PPPP_util::write(const Complex*,FILE*,bool))\n");
      exit(1);
    }
  }
  else{
    int length = fwrite(a, sizeof(Complex), 2*2 * 2*2, fp);
    if(length != 2*2 * 2*2){
      fprintf(stderr,
	      "ERROR(PPPP_util::write(const Complex*,FILE*,bool))\n");
      exit(1);
    }
  }
}

End_Name_Space_Yukawa
