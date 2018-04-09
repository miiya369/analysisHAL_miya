//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Potential
 * @brief   Definitions of the namespace for functions of potential
 * @author  Takaya Miyamoto
 * @since   Mon Nov  7 16:14:42 JST 2016
 */
//--------------------------------------------------------------------------

#ifndef POTENTIAL_HAL_H
#define POTENTIAL_HAL_H

#include "AnalysisHAL.h"
#include "ComplexField_Base.h"

namespace Potential {
   //! Return effective central potential (only 1st time derivative)
   ComplexField_BASE get_potential_T1(const ComplexField_BASE&,
                                      const ComplexField_BASE&,
                                      const double, const double);
   ComplexField_BASE get_potential_T1(const ComplexField_BASE&,
                                      const ComplexField_BASE&,
                                      const ComplexField_BASE&,
                                      const double, const double);
   //! Return effective central potential (with 2nd time derivative)
   ComplexField_BASE get_potential_T2(const ComplexField_BASE&,
                                      const ComplexField_BASE&,
                                      const ComplexField_BASE&,
                                      const double, const double);
   
   //! Return numerator of the eff-central potential (only 1st time derivative)
   ComplexField_BASE get_potential_T1_nume(const ComplexField_BASE&,
                                           const ComplexField_BASE&,
                                           const double, const double);
   ComplexField_BASE get_potential_T1_nume(const ComplexField_BASE&,
                                           const ComplexField_BASE&,
                                           const ComplexField_BASE&,
                                           const double, const double);
   //! Return numerator of the eff-central potential (with 2nd time derivative)
   ComplexField_BASE get_potential_T2_nume(const ComplexField_BASE&,
                                           const ComplexField_BASE&,
                                           const ComplexField_BASE&,
                                           const double, const double);
   
   //! Calculate the coupled channel potential (only 1st time derivative)
   void calc_CCpotential_T1(      ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const  double*, const double*,
                            const cdouble*, const int);
   void calc_CCpotential_T1(      ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const  double*, const double*,
                            const cdouble*, const int);
   //! Calculate the Coupled channel potential (with 2nd time derivative)
   void calc_CCpotential_T2(      ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const ComplexField_BASE*,
                            const  double*, const double*,
                            const cdouble*, const int);
   
   //! Return [S_{12} Psi]_{alpha,beta} for tensor force (Octet baryon)
   //! Require that the alpha == 4
   ComplexField_BASE get_S12_Psi   (const ComplexField_BASE&);
   //! Return [YD^star * Psi]_{alpha,beta} for tensor force (Octet baryon)
   //! Require that the alpha == 4
   ComplexField_BASE get_YDstar_Psi(const ComplexField_BASE&, const int);
   
   //! Calculate the central & tensor potential (only 1st time derivative)
   void calc_tensor_pot_T1(      ComplexField_BASE&, ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const double, const double);
   void calc_tensor_pot_T1(      ComplexField_BASE&, ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const double, const double);
   //! Calculate the central & tensor potential (with 2nd time derivative)
   void calc_tensor_pot_T2(      ComplexField_BASE&, ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const ComplexField_BASE&,
                           const double, const double);
   
   //! Calculate the central & tensor CC-potential (with 2nd time derivative)
   void calc_tensor_CCpot_T2(      ComplexField_BASE*, ComplexField_BASE*,
                             const ComplexField_BASE*,
                             const ComplexField_BASE*,
                             const ComplexField_BASE*,
                             const  double*, const double*,
                             const cdouble*, const int);
}

#endif
