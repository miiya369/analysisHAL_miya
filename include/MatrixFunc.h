//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Matrix
 * @brief   Header file of the function for matrix data
 * @author  Takaya Miyamoto
 * @since   Mon May 14 21:41:22 JST 2018
 */
//--------------------------------------------------------------------------

#ifndef MATRIX_FUNC_H
#define MATRIX_FUNC_H

//--------------------------------------------------------------------------
/**
 * @brief The namespace of the function for matrix data
 */
//--------------------------------------------------------------------------
namespace matrix_func
{
   //! I assume that X = double or complex<double>
   template <typename X>
   void solve_mat_gauss(const X*, X*, const X*, const int);
   
   template <typename X>
   void inverse_matrix (const X*, X*, const int);
}

#define ij(i,j) ((j) + Ndim * (i))

//--------------------------------------------------------------------------
/**
 * @brief The function for calculate AX=B using Gauss–Jordan elimination
 */
//--------------------------------------------------------------------------
template <typename X>
void matrix_func::solve_mat_gauss(const X *Amat, X *Xmat, 
				  const X *Bmat, const int Ndim) {
   X  tmp;
   X *mat = new X[Ndim*Ndim];
   
   for (int i=0; i<Ndim*Ndim; i++) {
      mat [i] = Amat[i];
      Xmat[i] = Bmat[i];
   }
   for (int i=0; i<Ndim; i++) { // A -> U
      tmp = mat[ij(i,i)];
      
      for (int j=0; j<Ndim; j++) {
	 Xmat[ij(i,j)] /= tmp;
	 mat [ij(i,j)] /= tmp;
      }
      for (int ii=i+1; ii<Ndim; ii++) {
	 tmp = mat[ij(ii,i)];
	 for (int jj=0; jj<Ndim; jj++) {
	    Xmat[ij(ii,jj)] -= Xmat[ij(i,jj)] * tmp;
	    mat [ij(ii,jj)] -= mat [ij(i,jj)] * tmp;
	 }
      }
   }
   for (   int  i=Ndim-1;  i!=0;  i--) // U -> E
      for (int ii=i;      ii!=0; ii--) {
	 tmp = mat[ij(ii-1,i)];
	 for (int jj=0; jj<Ndim; jj++) {
	    Xmat[ij(ii-1,jj)] -= Xmat[ij(i,jj)] * tmp;
	    mat [ij(ii-1,jj)] -= mat [ij(i,jj)] * tmp;
	 }
      }
   delete [] mat;
}

template <typename X>
void matrix_func::inverse_matrix(const X *imat, X *inv_mat, 
				 const int Ndim) {
   X *mat_E = new X[Ndim*Ndim];
   for (   int i=0; i<Ndim; i++) 
      for (int j=0; j<Ndim; j++) 
	 i==j ? mat_E[ij(i,j)] = 1.0 : mat_E[ij(i,j)] = 0.0;
   
   solve_mat_gauss(imat, inv_mat, mat_E, Ndim);
   
   delete [] mat_E;
}

#undef ij

#endif
