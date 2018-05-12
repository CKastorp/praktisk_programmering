#include"stdio.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"
void matrixtranspose(gsl_matrix* A, gsl_matrix* AT);

int QRsolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
int m=Q->size2;
int n=Q->size1;

if(b->size !=n) return GSL_EINVAL;
if(x->size!=m) return GSL_EINVAL;
if (R->size1!=R->size2)return GSL_ENOTSQR;
if (R->size1!=m)return GSL_EINVAL;
gsl_matrix* QT=gsl_matrix_alloc(m,n);
gsl_vector* RHS=gsl_vector_alloc(m);
matrixtranspose(Q,QT);

gsl_blas_dgemv(CblasNoTrans,1,QT,b,0,RHS);

int i,j;
for(i=m-1;i>=0;i--){
    double xi=gsl_vector_get(RHS,i)/gsl_matrix_get(R,i,i);
gsl_vector_set(x,i,xi);
for (j=0;j<i;j++){
    double RHS_new=gsl_vector_get(RHS,j)-gsl_matrix_get(R,j,i)*xi;
    gsl_vector_set(RHS,j,RHS_new);
}

}
gsl_vector_free(RHS);
gsl_matrix_free(QT);
    return GSL_SUCCESS;
}