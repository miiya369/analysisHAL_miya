//--------------------------------------------------------------------------
/**
 * @file
 * @brief   representations of cubic groups and related topics
 * @ingroup Cubic_Group
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Wed Aug 18 19:30:41 JST 2010
 *
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_CHARACTER_H
#define IS_INCLUDED_CHARACTER_H

#include <complex>
#include <yukawa/yukawa.h>
#include <yukawa/matrix.h>
#include <yukawa/GPM.h>

Start_Name_Space_Yukawa

namespace cubic_group {
  enum Representation_O  {A1=0, A2=1, E=2, T2=3, T1=4, UNKNOWN=-1};
  enum Representation_C4 {A=0,  B=1,  E1=2, E2=3};
  enum Parity            {Plus=0, Minus=1};
  enum Spin              {Singlet=0, Triplet=1};

  extern std::complex<double> character_O[];
  extern std::complex<double> character_C4[];

  inline
  std::complex<double>& chara_O(int i, int rep) {return character_O[i+5*rep];}

  inline
  std::complex<double>& chara_C4(int i,int rep) {return character_C4[i+4*rep];}
  
  extern matrix<std::complex<double> >
  projection_matrix_O(const std::complex<double>          chara[],
		      const matrix<std::complex<double> > mat[]);

  extern matrix<std::complex<double> >
  projection_matrix_C4(const std::complex<double>          chara[],
		       const matrix<std::complex<double> > mat[]);
}

End_Name_Space_Yukawa

#endif
