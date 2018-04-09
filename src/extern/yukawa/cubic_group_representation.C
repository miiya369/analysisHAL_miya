//--------------------------------------------------------------------------
/**
 * @file
 * @brief   representations of cubic groups and related topics
 * @ingroup Cubic_Group
 * @author  N.Ishii
 * @since   Wed Aug 18 19:36:15 JST 2010
 */
//--------------------------------------------------------------------------

#include <complex>
#include <yukawa/yukawa.h>
#include <yukawa/cubic_group_representation.h>
#include <yukawa/matrix.h>

typedef std::complex<double> Complex;

Start_Name_Space_Yukawa

typedef matrix<Complex> Matrix;

//--------------------------------------------------------------------------
/**
 * @brief character of the cubic group representations
 *
 * Complex cubic_group::character_O[]={  
 *   //E,  8C3,  3C2,  6C2,  6C4       
 *   1.0,  1.0,  1.0,  1.0,  1.0, // A1
 *   1.0,  1.0,  1.0, -1.0, -1.0, // A2
 *   2.0, -1.0,  2.0,  0.0,  0.0, // E 
 *   3.0,  0.0, -1.0,  1.0, -1.0, // T2
 *   3.0,  0.0, -1.0, -1.0,  1.0  // T1
 * };                                  
 *
 */
//--------------------------------------------------------------------------

Complex cubic_group::character_O[]={
  //E,  8C3,  3C2,  6C2,  6C4
  1.0,  1.0,  1.0,  1.0,  1.0, // A1
  1.0,  1.0,  1.0, -1.0, -1.0, // A2
  2.0, -1.0,  2.0,  0.0,  0.0, // E
  3.0,  0.0, -1.0,  1.0, -1.0, // T2
  3.0,  0.0, -1.0, -1.0,  1.0  // T1
};                                  

//--------------------------------------------------------------------------
/**
 * @brief characters of C4 group
 *
 * Complex cubic_group::character_C4[]={  
 *   //E   C4    C2    C4^3               
 *   1.0,  1.0,  1.0,  1.0, // A          
 *   1.0, -1.0,  1.0, -1.0, // B          
 *   1.0,  I,   -1.0, -I,   // E+         
 *   1.0, -I,   -1.0,  I    // E-         
 * };                                     
 *
 */
//--------------------------------------------------------------------------

Complex cubic_group::character_C4[]={
  //E   C4                   C2    C4^3
  1.0,  1.0,                 1.0,  1.0,                // A
  1.0, -1.0,                 1.0, -1.0,                // B
  1.0,  Complex(0.0,1.0),   -1.0, -Complex(0.0,1.0),   // E+
  1.0, -Complex(0.0,1.0),   -1.0,  Complex(0.0,1.0)    // E-
};                                     

//--------------------------------------------------------------------------
/**
 * @brief projection matrix for O
 */
//--------------------------------------------------------------------------

Matrix
cubic_group::projection_matrix_O(const Complex chara[], const Matrix mat[])
{
  Matrix ret(mat[0].size());

  // E
  for(int ith = 0; ith < 1; ith++){
    ret += conj(chara[0]) * mat[ith];
  }

  // 6C4
  for(int ith = 1; ith < 7; ith++){
    ret += conj(chara[4]) * mat[ith];
  }

  // 8C3
  for(int ith = 7; ith < 15; ith++){
    ret += conj(chara[1]) * mat[ith];
  }

  // 3C2
  for(int ith = 15; ith < 18; ith++){
    ret += conj(chara[2]) * mat[ith];
  }

  // 6C2
  for(int ith = 18; ith < 24; ith++){
    ret += conj(chara[3]) * mat[ith];
  }

  return ret * real(chara[0])/24.0;  
}

//--------------------------------------------------------------------------
/**
 * @brief projection matrix for C4
 */
//--------------------------------------------------------------------------

Matrix
cubic_group::projection_matrix_C4(const Complex chara[], const Matrix mat[])
{
  Matrix ret(mat[0].size());

  // E
  ret += conj(chara[0]) * mat[0];

  // C4 (c_z)
  ret += conj(chara[1]) * mat[1];

  // C2 (c_z^2)
  ret += conj(chara[2]) * mat[2];

  // C4^3 (c_z^3)
  ret += conj(chara[3]) * mat[3];

  return ret * real(chara[0])/4.0;
}

End_Name_Space_Yukawa

