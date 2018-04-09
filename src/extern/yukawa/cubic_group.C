//--------------------------------------------------------------------------
/**
 * @file    
 * @brief   Cubic rotation group and its (Pauli) representation matrix
 * @ingroup Cubic_Group
 * @author  N.Ishii
 * @date    Fri Aug 24 16:04:27 JST 2007
 */
//--------------------------------------------------------------------------

#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/cubic_group.h>
#include <yukawa/GPM.h>
#include <yukawa/matrix.h>

Start_Name_Space_Yukawa

typedef std::complex<double> Complex;
typedef matrix<Complex>      Matrix;

/// cubic transformation group
GPM<int> cubic_group::cube[24];

/// Pauli matrix
matrix<Complex> cubic_group::S[24];

static int is_initialized=0;

//--------------------------------------------------------------------------
/**
 * Initialization of the 24 elements of the cubic transformation group
 */
//--------------------------------------------------------------------------
void cubic_group::initialize()
{
  if(is_initialized==1) return;

  using namespace cubic_group;

  //---------------------------------------------------------------------
  // cubic rotation group (orbit)
  //---------------------------------------------------------------------

  typedef GPM<int> GPMI;

  //GPMI& e( cube[0]);
  //GPMI& Cx(cube[1]);
  //GPMI& Cy(cube[2]);
  //GPMI& Cz(cube[3]);

  cube[0]  = GPMI(0,1,2, 1,1,1);
  //--------------------------------------------------------------------------
  cube[2]  = GPMI(2,1,0, 1,1,-1);
  cube[3]  = GPMI(1,0,2, -1,1,1);

  construct_representation_matrices(cube);

  //---------------------------------------------------------------------
  // cubic rotation group (Pauli spinor)
  //---------------------------------------------------------------------

  const double sqrt2 = 1.41421356237309504880;

  for(int i = 0; i < 24; i++) S[i] = Matrix(2);

  S[0](0,0) = 1.0;  S[0](0,1) = 0.0;
  S[0](1,0) = 0.0;  S[0](1,1) = 1.0;

  S[2](0,0) = 1.0/sqrt2; S[2](0,1) = -1.0/sqrt2;
  S[2](1,0) = 1.0/sqrt2; S[2](1,1) =  1.0/sqrt2;

  S[3](0,0) = Complex(1.0,-1.0)/sqrt2;  S[3](0,1) = 0.0;
  S[3](1,0) = 0.0;                      S[3](1,1) = Complex(1.0,1.0)/sqrt2;
  
  construct_representation_matrices(S);

  is_initialized=1;
}

End_Name_Space_Yukawa
