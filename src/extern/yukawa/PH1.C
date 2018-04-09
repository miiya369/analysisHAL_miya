//--------------------------------------------------------------------------
/**
 * @file
 * @brief
 * @ingroup PH1_group
 * @author  N.Ishii
 * @since   Fri May 21 21:36:04 JST 2010
 */
//--------------------------------------------------------------------------

#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/PH1.h>
#include <yukawa/GPM.h>
#include <yukawa/matrix.h>
#include <yukawa/cubic_group.h>
#include <yukawa/misc.h>

#if TARGET == YUKAWA_BGQ
#define BGQ_130720
#endif

#ifdef BGQ_130720
#include <builtins.h>
#define Flip(R)                vec_sldw(R,R,2)
#define VEC_xmadd(   v1,v2,v3) vec_xmadd(   v1,v2,v3)
#define VEC_xxnpmadd(v1,v2,v3) vec_xxnpmadd(v2,v1,v3)
#endif

Start_Name_Space_Yukawa

typedef std::complex<double> Complex;
typedef matrix<Complex>      Matrix;
typedef GPM<double>          GPMR;

static bool is_initialized = false;

static GPMR CxCxCxC;

static Matrix SxSxSstarxSstar[24];

static GPMR UxUxUxU;

void PH1::initialize()
{
  if (is_initialized == true) return;

  //--------------------------------------------------------------------------
  // charge conjugation matrix
  //--------------------------------------------------------------------------

  GPMR C(1,0, -1.0,1.0);
  GPMR CxC     = Kronecker_product(C,C);
  
  CxCxCxC = Kronecker_product(CxC,CxC);

  //--------------------------------------------------------------------------
  // spinor roration
  //--------------------------------------------------------------------------

  Matrix S[24];
  for(int i = 0; i < 24; i++) S[i] = Matrix(2);
  const double sqrt2 = 1.41421356237309504880;

  S[0](0,0) = 1.0; S[0](0,1) = 0.0;
  S[0](1,0) = 0.0; S[0](1,1) = 1.0;
  
  S[2](0,0) = 1.0/sqrt2; S[2](0,1) = -1.0/sqrt2;
  S[2](1,0) = 1.0/sqrt2; S[2](1,1) =  1.0/sqrt2;
  
  S[3](0,0) = Complex(1.0,-1.0)/sqrt2; S[3](0,1) = 0.0;
  S[3](1,0) = 0.0;                     S[3](1,1) = Complex(1.0,1.0)/sqrt2;

  cubic_group::construct_representation_matrices(S);

  Matrix Sstar[24];
  for(int i = 0; i < 24; i++) Sstar[i] = Matrix(2);
  for(int i = 0; i < 24; i++){
    Sstar[i](0,0) = conj(S[i](0,0)); Sstar[i](0,1) = conj(S[i](0,1));
    Sstar[i](1,0) = conj(S[i](1,0)); Sstar[i](1,1) = conj(S[i](1,1));
  }

  Matrix SxS[24];
  Matrix SstarxSstar[24];
  
  for(int i = 0; i < 24; i++){
    SxS[i]         = Kronecker_product(S[i],S[i]);
    SstarxSstar[i] = Kronecker_product(Sstar[i],Sstar[i]);
  }
  
  for(int i = 0; i < 24; i++)
    SxSxSstarxSstar[i] = Kronecker_product(SxS[i], SstarxSstar[i]);

#if 0
  FILE* fp;
  if((fp = fopen("PH1.SxSxSstarxSstar.txt","w"))==NULL){
    fprintf(stderr, "cannot fopen 'PH1.SxSxSstarxSstar.txt'\n");
    exit(1);
  }
  for(int ith = 0; ith < 24; ith++){
    Matrix mat(SxSxSstarxSstar[ith]);
    int dim = mat.size();
    fprintf(fp, "### ith=%d\n", ith);
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++)
	fprintf(fp, "(%1.2f, %1.2f)  ", real(mat(i,j)), imag(mat(i,j)));
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
#endif

  //--------------------------------------------------------------------------
  // cps gamma matrix convention
  //--------------------------------------------------------------------------

  GPMR U(1,0, 1.0,-1.0);
  GPMR UxU = Kronecker_product(U,U);

  UxUxUxU  = Kronecker_product(UxU,UxU);

  is_initialized = true;
}

//--------------------------------------------------------------------------
/**
 *
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

PH1& PH1::charge_conjugation()
{
  Complex tmp[2*2 * 2*2];
  for(int i = 0; i < 16; i++)
    tmp[i] = CxCxCxC.Z[i] * conj(a[CxCxCxC.I[i]]);

  for(int i = 0; i < 2 * 16; i++)
    ((double*)a)[i] = ((double*)tmp)[i];

  return *this;
}

//---------------------------------------------------------------------
/**
 * @brief spinor rotation
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
//---------------------------------------------------------------------

PH1& PH1::rotate(int i_th)
#ifdef BGQ_130720
{
#define Buf(n) buf_ ## n ## _
#define Mat(n) mat_ ## n ## _
#define Ret(n) ret_ ## n ## _

  Complex (*mat)[16] = (Complex (*)[16])&SxSxSstarxSstar[i_th](0,0);

#pragma disjoint(*a,*mat)

  __alignx(32,a);
  __alignx(32,mat);

  vector4double Buf(0), Buf(2),  Buf(4),  Buf(6);
  vector4double Buf(8), Buf(10), Buf(12), Buf(14);

  vector4double Mat(0), Mat(2),  Mat(4),  Mat(6);
  vector4double Mat(8), Mat(10), Mat(12), Mat(14);

  vector4double Ret(0), Ret(2),  Ret(4),  Ret(6);
  vector4double Ret(8), Ret(10), Ret(12), Ret(14);

  vector4double A, B;

  vector4double Zero=(vector4double)(0.0, 0.0, 0.0, 0.0);
  vector4double One =(vector4double)(1.0, 1.0, 1.0, 1.0);

  Buf(0)  = vec_ld(0, (double*)&((Complex*)a)[0]);
  Buf(2)  = vec_ld(0, (double*)&((Complex*)a)[2]);
  Buf(4)  = vec_ld(0, (double*)&((Complex*)a)[4]);
  Buf(6)  = vec_ld(0, (double*)&((Complex*)a)[6]);
  Buf(8)  = vec_ld(0, (double*)&((Complex*)a)[8]);
  Buf(10) = vec_ld(0, (double*)&((Complex*)a)[10]);
  Buf(12) = vec_ld(0, (double*)&((Complex*)a)[12]);
  Buf(14) = vec_ld(0, (double*)&((Complex*)a)[14]);

  for(int i = 0; i < 16; i+=2){
    //---------------------------------------------------------------------
    // zenhan
    //---------------------------------------------------------------------
    Mat(0)  = vec_ld(0, (double*)&mat[i][0]);
    Mat(2)  = vec_ld(0, (double*)&mat[i][2]);
    Mat(4)  = vec_ld(0, (double*)&mat[i][4]);
    Mat(6)  = vec_ld(0, (double*)&mat[i][6]);
    Mat(8)  = vec_ld(0, (double*)&mat[i][8]);
    Mat(10) = vec_ld(0, (double*)&mat[i][10]);
    Mat(12) = vec_ld(0, (double*)&mat[i][12]);
    Mat(14) = vec_ld(0, (double*)&mat[i][14]);


    Ret(0)  = VEC_xmadd(   Mat(0),  Buf(0),  Zero);
    Ret(2)  = VEC_xmadd(   Mat(2),  Buf(2),  Zero);
    Ret(4)  = VEC_xmadd(   Mat(4),  Buf(4),  Zero);
    Ret(6)  = VEC_xmadd(   Mat(6),  Buf(6),  Zero);
    Ret(8)  = VEC_xmadd(   Mat(8),  Buf(8),  Zero);
    Ret(10) = VEC_xmadd(   Mat(10), Buf(10), Zero);
    Ret(12) = VEC_xmadd(   Mat(12), Buf(12), Zero);
    Ret(14) = VEC_xmadd(   Mat(14), Buf(14), Zero);

    Ret(0)  = VEC_xxnpmadd(Mat(0),  Buf(0),  Ret(0) );
    Ret(2)  = VEC_xxnpmadd(Mat(2),  Buf(2),  Ret(2) );
    Ret(4)  = VEC_xxnpmadd(Mat(4),  Buf(4),  Ret(4) );
    Ret(6)  = VEC_xxnpmadd(Mat(6),  Buf(6),  Ret(6) );
    Ret(8)  = VEC_xxnpmadd(Mat(8),  Buf(8),  Ret(8) );
    Ret(10) = VEC_xxnpmadd(Mat(10), Buf(10), Ret(10));
    Ret(12) = VEC_xxnpmadd(Mat(12), Buf(12), Ret(12));
    Ret(14) = VEC_xxnpmadd(Mat(14), Buf(14), Ret(14));


    Ret(0)  = vec_madd(One, Ret(0),  Ret(2));
    Ret(4)  = vec_madd(One, Ret(4),  Ret(6));
    Ret(8)  = vec_madd(One, Ret(8),  Ret(10));
    Ret(12) = vec_madd(One, Ret(12), Ret(14));

    Ret(0)  = vec_madd(One, Ret(0),  Ret(4));
    Ret(8)  = vec_madd(One, Ret(8),  Ret(12));

    Ret(0)  = vec_madd(One, Ret(0),  Ret(8));
    Ret(2)  = Flip(Ret(0));
    A       = vec_madd(One, Ret(0),  Ret(2));

    //---------------------------------------------------------------------
    // kohan
    //---------------------------------------------------------------------

    Mat(0)  = vec_ld(0, (double*)&mat[i+1][0]);
    Mat(2)  = vec_ld(0, (double*)&mat[i+1][2]);
    Mat(4)  = vec_ld(0, (double*)&mat[i+1][4]);
    Mat(6)  = vec_ld(0, (double*)&mat[i+1][6]);
    Mat(8)  = vec_ld(0, (double*)&mat[i+1][8]);
    Mat(10) = vec_ld(0, (double*)&mat[i+1][10]);
    Mat(12) = vec_ld(0, (double*)&mat[i+1][12]);
    Mat(14) = vec_ld(0, (double*)&mat[i+1][14]);


    Ret(0)  = VEC_xmadd(   Mat(0),  Buf(0),  Zero);
    Ret(2)  = VEC_xmadd(   Mat(2),  Buf(2),  Zero);
    Ret(4)  = VEC_xmadd(   Mat(4),  Buf(4),  Zero);
    Ret(6)  = VEC_xmadd(   Mat(6),  Buf(6),  Zero);
    Ret(8)  = VEC_xmadd(   Mat(8),  Buf(8),  Zero);
    Ret(10) = VEC_xmadd(   Mat(10), Buf(10), Zero);
    Ret(12) = VEC_xmadd(   Mat(12), Buf(12), Zero);
    Ret(14) = VEC_xmadd(   Mat(14), Buf(14), Zero);

    Ret(0)  = VEC_xxnpmadd(Mat(0),  Buf(0),  Ret(0) );
    Ret(2)  = VEC_xxnpmadd(Mat(2),  Buf(2),  Ret(2) );
    Ret(4)  = VEC_xxnpmadd(Mat(4),  Buf(4),  Ret(4) );
    Ret(6)  = VEC_xxnpmadd(Mat(6),  Buf(6),  Ret(6) );
    Ret(8)  = VEC_xxnpmadd(Mat(8),  Buf(8),  Ret(8) );
    Ret(10) = VEC_xxnpmadd(Mat(10), Buf(10), Ret(10));
    Ret(12) = VEC_xxnpmadd(Mat(12), Buf(12), Ret(12));
    Ret(14) = VEC_xxnpmadd(Mat(14), Buf(14), Ret(14));


    Ret(0)  = vec_madd(One, Ret(0),  Ret(2));
    Ret(4)  = vec_madd(One, Ret(4),  Ret(6));
    Ret(8)  = vec_madd(One, Ret(8),  Ret(10));
    Ret(12) = vec_madd(One, Ret(12), Ret(14));

    Ret(0)  = vec_madd(One, Ret(0),  Ret(4));
    Ret(8)  = vec_madd(One, Ret(8),  Ret(12));

    Ret(0)  = vec_madd(One, Ret(0),  Ret(8));
    Ret(2)  = Flip(Ret(0));
    B       = vec_madd(One, Ret(0),  Ret(2));

    //---------------------------------------------------------------------
    // atoshori
    //---------------------------------------------------------------------

    A = vec_sldw(A,B,2);

    vec_st(A, 0, (double*)&((Complex*)a)[i]);

  }

  return *this;

#undef Ret
#undef Mat
#undef Buf
}
#else
{
  Matrix& mat(SxSxSstarxSstar[i_th]);

  Complex tmp[16];
  for(int i = 0; i < 2 * 16; i++)
    ((double*)tmp)[i] = ((double*)a)[i];

  for(int i = 0; i < 16; i++){
    Complex sum = 0.0;
    for(int j = 0; j < 16; j++) sum += mat(i,j) * tmp[j];
    a[i] = sum;
  }
  return *this;
}
#endif

/// read
PH1& PH1::read(FILE* fp, bool is_change_endian)
{
  if(fread(a, sizeof(Complex), 2*2 * 2*2, fp) != 16){
    fprintf(stderr, "ERROR(PH1::read)\n");
    exit(1);
  }
  if(is_change_endian == true){
    yukawa::change_endian((double*)a, 2 * 2*2 * 2*2);
  }
  return *this;
}

/// write
const
PH1& PH1::write(FILE* fp, bool is_change_endian) const
{
  if(is_change_endian == true){
    Complex b[2*2 * 2*2];
    for(int i = 0; i < 2*2 * 2*2; i++) b[i] = a[i];
    yukawa::change_endian((double*)b, 2 * 2*2 * 2*2);
    if(fwrite(b, sizeof(Complex), 2*2 * 2*2, fp) != 16){
      fprintf(stderr, "ERROR(PH1::write)\n");
      exit(1);
    }
  }
  else{
    if(fwrite(a, sizeof(Complex), 2*2 * 2*2, fp) != 16){
      fprintf(stderr, "ERROR(PH1::write)\n");
      exit(1);
    }
  }
  return *this;
}

/// endian change
PH1& PH1::change_endian()
{
  char tmp[8];
  for(int i = 0; i < 2 * 2*2 * 2*2; i++){
    double* ptr = &((double*)a)[i];
    for(int j = 0; j < 8; j++) tmp[j] = ((char*)ptr)[j];
    for(int j = 0; j < 8; j++) ((char*)ptr)[j] = tmp[7 - j];
  }
  return *this;
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

PH1& PH1::cps2standard()
{
  Complex tmp[16];
  Complex* ptr = (Complex*)&(*this)(0,0,0,0);

  for(int i = 0; i < 16; i++)
    tmp[i] = UxUxUxU.Z[i] * ptr[UxUxUxU.I[i]];

    for(int i = 0; i < 2 * 16; i++)
      ((double*)ptr)[i] = ((double*)tmp)[i];
    return *this;
}

/// helper function version
PH1& PH1_cps2standard(PH1& obj)
{
  Complex tmp[16];
  Complex* ptr = &obj(0,0,0,0);

  for(int i = 0; i < 16; i++)
    tmp[i] = UxUxUxU.Z[i] * ptr[UxUxUxU.I[i]];

    for(int i = 0; i < 2 * 16; i++)
      ((double*)ptr)[i] = ((double*)tmp)[i];
    return obj;
}

End_Name_Space_Yukawa

