//--------------------------------------------------------------------------
/**
 * @file
 *
 * @brief   Cubic transformation group
 * @ingroup Cubic_Group
 * @ingroup Header_File
 * @author  N.Ishii
 * @date    Sat Feb 16 21:58:31 JST 2008
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup Cubic_Group cubic group related topics
 * @ingroup  General_Util
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_CUBIC_GROUP_H
#define IS_INCLUDED_CUBIC_GROUP_H

#include <yukawa/yukawa.h>
#include <yukawa/GPM.h>
#include <yukawa/matrix.h>

#include <yukawa/cubic_group_representation.h>

Start_Name_Space_Yukawa

namespace cubic_group {

  //---------------------------------------------------------------------
  /**
   * @brief template function to construct an array of type @a T consisting
   *        of 24 elements of representation matrices
   *        for the cubic transformation group.
   *
   * It is assumed that type @T has operations<br>
   * (i)  multiplication: T operator*(T lhs,const T& rhs)<br>
   * (ii) inverse: T dagger(const T& obj)<br>
   *
   * Provide the identity element  and the elements corresponding to the
   * rotation  by 90  degree around  y, z  axis through  ary[0], ary[2],
   * ary[3], respectively. The element  corresponding to the rotation of
   * 90 degree  around x axis is  generated in ary[1], and  the other 20
   * elements are generated in ary[n] (n=4,5, ..., 23).
   *
   *
   */
  //---------------------------------------------------------------------
  
  template <class T> inline
    void construct_representation_matrices(T* ary)
    {
      T& Cx(ary[1]);
      T& Cy(ary[2]);
      T& Cz(ary[3]);
    
      T& CxInv(ary[4]);
      T& CyInv(ary[5]);
      T& CzInv(ary[6]);

      CyInv  = Cy * Cy * Cy; //dagger(Cy);
      CzInv  = Cz * Cz * Cz; // dagger(Cz);

      Cx     = Cy * Cz * dagger(Cy);
      CxInv  = Cx * Cx * Cx; //dagger(Cx);
      //-----------------------------------------------------------------------
      ary[7]  = Cx * Cy;
      ary[8]  = Cx * Cy  *  Cx * Cy;
    
      ary[9]  = Cy * Cx;
      ary[10] = Cz * Cy;
      ary[11] = Cx * Cz;
    
      ary[12] = Cy * Cx  *  Cy * Cx;
      ary[13] = Cz * Cy  *  Cz * Cy;
      ary[14] = Cx * Cz  *  Cx * Cz;
      //-----------------------------------------------------------------------
      ary[15] = Cx * Cx;
      ary[16] = Cy * Cy;
      ary[17] = Cz * Cz;
    
      ary[18] = Cx * Cy * Cz;
      ary[19] = Cy * Cz * Cx;
      ary[20] = Cz * Cx * Cy;
    
      ary[21] = CxInv * CyInv * Cz;
      ary[22] = CyInv * CzInv * Cx;
      ary[23] = CzInv * CxInv * Cy;
    }

  //--------------------------------------------------------------------------
  /**
   * @brief C4 subgroup
   *
   *
   */
  //--------------------------------------------------------------------------
  template <class T> inline
  void pickup_C4_subgroup(const T* ary, T* ret)
  {
    ret[0] = ary[0];  // E
    ret[1] = ary[3];  // C_z
    ret[2] = ary[17]; // C_z^2
    ret[3] = ary[6];  // C_z^3
  }

  extern GPM<int>                      cube[24];
  extern matrix<std::complex<double> > S[24];

  extern void initialize();
}

End_Name_Space_Yukawa

#endif
