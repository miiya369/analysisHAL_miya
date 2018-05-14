#include <CommonIncl.h>
#include <MatrixFunc.h>

int main() {
  cdouble Amat[2*2] = {1, -4,
		       1,  2};
  
  cdouble Xmat[2*2];
  
  for (int i=0; i<2*2; i++) Xmat[i] = Amat[i];
  
  for (  int i2=0; i2<2; i2++) {
    for (int j2=0; j2<2; j2++) printf("%8.4lf ", Xmat[j2+2*i2].real());
    printf("\n");
  } printf("\n");
  
  matrix_func::inverse_matrix(Amat, Xmat, 2);
  
  for (  int i2=0; i2<2; i2++) {
    for (int j2=0; j2<2; j2++) printf("%8.4lf ", Xmat[j2+2*i2].real());
    printf("\n");
  } printf("\n");
  
  cdouble Bmat[2*2];
  for (  int i=0; i<2; i++)
    for (int j=0; j<2; j++) {
      Bmat[j+2*i] = 0.0;
      
      for (int k=0; k<2; k++)
	Bmat[j+2*i] += Xmat[k+2*i] * Amat[j+2*k];
    }
  
  for (  int i2=0; i2<2; i2++) {
    for (int j2=0; j2<2; j2++) printf("%8.4lf ", Bmat[j2+2*i2].real());
    printf("\n");
  } printf("\n");
  
  return 0;
}
