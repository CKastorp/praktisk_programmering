#include"gsl/gsl_errno.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"math.h"
int QRfactor(gsl_matrix*Q,gsl_matrix*R);
int QRsolve(gsl_matrix*Q,gsl_matrix*R,gsl_vector*b,gsl_vector*x);

int QRleastsquare(gsl_vector*x,gsl_vector*y,gsl_vector*dy,gsl_vector*c,double fitfunction(double x,int i)){
int n=x->size,k=c->size;
if(y->size!=n)return GSL_EINVAL;
if(dy->size!=n)return GSL_EINVAL;

gsl_matrix* A=gsl_matrix_alloc(n,k);
gsl_matrix* R=gsl_matrix_alloc(k,k);
for (int j=0;j<k;j++){
for (int i=0;i<n;i++){double xi=gsl_vector_get(x,i);
 gsl_matrix_set(A,i,j,fitfunction(xi,j));}
}
QRfactor(A,R);
QRsolve(A,R,y,c);

gsl_matrix_free(A);
gsl_matrix_free(R);
    return GSL_SUCCESS;
}