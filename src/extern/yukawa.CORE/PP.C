//--------------------------------------------------------------------------
/**
 * @file
 * @brief
 * @ingroup PP_group
 * @author  N.Ishii
 * @since   Sun Aug 15 19:55:59 JST 2010
 */
//--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/4pt.h>
#include <yukawa/matrix.h>
#include <yukawa/cubic_group.h>
#include <yukawa.CORE/PP.h>
#include <yukawa.CORE/PPPP.h>
#include <yukawa/misc.h>



Start_Name_Space_Yukawa

//inline double conj(double x) { return x; }

typedef std::complex<double> Complex;
typedef matrix<Complex>      Matrix;

Matrix PP::S[24];
Matrix PP::SxS[24];

GPM<double> PP::C;
GPM<double> PP::CxC;

double PP::character_O_ary[] =
  {1.0,  1.0,  1.0,  1.0,  1.0,
   1.0,  1.0,  1.0, -1.0, -1.0,
   2.0, -1.0,  2.0,  0.0,  0.0,
   3.0,  0.0, -1.0,  1.0, -1.0,
   3.0,  0.0, -1.0, -1.0,  1.0};

static bool is_initialized = false;

void PP::initialize()
{
  if(is_initialized == true) return;

  cubic_group::initialize();
  PPPP_util::initialize();

  C = PPPP_util::C;
  CxC = Kronecker_product(C,C);

  for(int i = 0; i < 24; i++){
    S[i] = PPPP_util::S[i];
  }

  for(int i = 0; i < 24; i++){
    SxS[i] = Kronecker_product(S[i], S[i]);
  }

  is_initialized = true;
}

PP& PP::rotate(int ith)
{
  Matrix& mat(SxS[ith]);

  Complex tmp[4];
  for(int i = 0; i < 2 * 4; i++)
    ((double*)tmp)[i] = ((double*)a)[i];

  for(int i = 0; i < 4; i++){
    Complex sum = 0.0;
    for(int j = 0; j < 4; j++) sum += mat(i,j) * tmp[j];
    a[i] = sum;
  }
  return *this;
}

PP& PP::spin_singlet()
{
  (*this)(0,0) = 0.0;
  (*this)(1,1) = 0.0;

  Complex val = 0.5*((*this)(0,1) - (*this)(1,0));

  (*this)(0,1) =  val;
  (*this)(1,0) = -val;

  return *this;
}

PP& PP::spin_triplet()
{
  Complex val = 0.5*((*this)(0,1) + (*this)(1,0));

  (*this)(0,1) = val;
  (*this)(1,0) = val;

  return *this;
}

//--------------------------------------------------------------------------
// additionals
//--------------------------------------------------------------------------

#if 0

//--------------------------------------------------------------------------
/**
 * @brief L projection
 */
//--------------------------------------------------------------------------

four_point<PP> L_projection(const four_point<PP>& four, int rep)
{
  typedef four_point<PP> Four_Point;

  if(rep < 0 || 5 <= rep){
    fprintf(stderr,
	    "ERROR(L_projection): invalid value of rep=%d\n", rep);
    exit(1);
  }

  const double *chi = &PP::character_O(0, rep);

  Four_Point sum(conj(chi[0]) * four);

  // 6C4
  for(int ith = 1; ith < 7; ith++){
    sum += conj(chi[4]) * rotate(four, ith, true);
  }

  // 8C3
  for(int ith = 7; ith < 15; ith++){
    sum += conj(chi[1]) * rotate(four, ith, true);
  }

  // 3C2
  for(int ith = 15; ith < 18; ith++){
    sum += conj(chi[2]) * rotate(four, ith, true);
  }

  // 6C2
  for(int ith = 18; ith < 24; ith++){
    sum += conj(chi[3]) * rotate(four, ith, true);
  }

  return (chi[0]/24.0) * sum;
}

//--------------------------------------------------------------------------
/**
 * @brief L projection (Nanpa)
 * @param crep = "A1", "A2", "E", "T2", "T1"
 */
//--------------------------------------------------------------------------

four_point<PP> L_projection(const four_point<PP>& four, const char* crep)
{
  int rep=(int)PP::UNKNOWN;

  if(strcmp(crep,"A1")==0){
    rep = (int)PP::A1;
  }
  else if(strcmp(crep,"A2")==0){
    rep = (int)PP::A2;
  }
  else if(strcmp(crep,"E")==0){
    rep = (int)PP::E;
  }
  else if(strcmp(crep,"T2")==0){
    rep = (int)PP::T2;
  }
  else if(strcmp(crep,"T1")==0){
    rep = (int)PP::T1;
  }
  else{
    fprintf(stderr,
	    "ERROR(L_projection): invalid crep = \"%s\"\n", crep);
    exit(1);
  }

  return L_projection(four, rep);
}

#endif

End_Name_Space_Yukawa
