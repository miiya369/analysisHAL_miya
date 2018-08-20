//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Potential
 * @brief   Implementations for single-channel effective central potentials
 * @author  Takaya Miyamoto
 * @since   Tue Feb  2 18:57:38 JST 2016
 */
//--------------------------------------------------------------------------

#include <Potential.h>

ComplexField_BASE Potential::get_potential_T1(const ComplexField_BASE &Rm,
                                              const ComplexField_BASE &R0,
                                              const double m1, const double m2) {
   DEBUG_LOG
   double mu = (m1 * m2) / (m1 + m2);
   
   return ( R0.lap()/(2.0*mu) /* Laplacian */ +
           (Rm - R0)          /* T1 deriv. */) / R0;
}
ComplexField_BASE Potential::get_potential_T1(const ComplexField_BASE &Rm,
                                              const ComplexField_BASE &R0,
                                              const ComplexField_BASE &Rp,
                                              const double m1, const double m2) {
   DEBUG_LOG
   double mu = (m1 * m2) / (m1 + m2);
   
   return ( R0.lap()/(2.0*mu) /* Laplacian */ +
           (Rm - Rp)/ 2.0     /* T1 deriv. */) / R0;
}
ComplexField_BASE Potential::get_potential_T2(const ComplexField_BASE &Rm,
                                              const ComplexField_BASE &R0,
                                              const ComplexField_BASE &Rp,
                                              const double m1, const double m2) {
   DEBUG_LOG
   double mu = (m1 * m2) / (m1 + m2);
   double dl = (m1 - m2) / (m1 + m2);
   
   return ( R0.lap()/(2.0*mu)                          /* Laplacian */ +
           (Rm - Rp)/ 2.0                              /* T1 deriv. */ +
           (Rm + Rp - 2.0*R0)*(1.0+3.0*dl*dl)/(8.0*mu) /* T2 deriv. */
           ) / R0;
}

ComplexField_BASE Potential::get_potential_T1_nume(const ComplexField_BASE &Rm,
                                                   const ComplexField_BASE &R0,
                                                   const double m1,
                                                   const double m2) {
   DEBUG_LOG
   double mu = (m1 * m2) / (m1 + m2);
   
   return ( R0.lap()/(2.0*mu) /* Laplacian */ +
           (Rm - R0)          /* T1 deriv. */);
}
ComplexField_BASE Potential::get_potential_T1_nume(const ComplexField_BASE &Rm,
                                                   const ComplexField_BASE &R0,
                                                   const ComplexField_BASE &Rp,
                                                   const double m1,
                                                   const double m2) {
   DEBUG_LOG
   double mu = (m1 * m2) / (m1 + m2);
   
   return ( R0.lap()/(2.0*mu) /* Laplacian */ +
           (Rm - Rp)/ 2.0     /* T1 deriv. */);
}
ComplexField_BASE Potential::get_potential_T2_nume(const ComplexField_BASE &Rm,
                                                   const ComplexField_BASE &R0,
                                                   const ComplexField_BASE &Rp,
                                                   const double m1,
                                                   const double m2) {
   DEBUG_LOG
   double mu = (m1 * m2) / (m1 + m2);
   double dl = (m1 - m2) / (m1 + m2);
   
   return ( R0.lap()/(2.0*mu)                          /* Laplacian */ +
           (Rm - Rp)/ 2.0                              /* T1 deriv. */ +
           (Rm + Rp - 2.0*R0)*(1.0+3.0*dl*dl)/(8.0*mu) /* T2 deriv. */
           );
}
