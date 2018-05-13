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

gsl_vector* RHS=gsl_vector_alloc(m);

gsl_blas_dgemv(CblasTrans,1,Q,b,0,RHS);

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

    return GSL_SUCCESS;
}

int QRinverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){

int n=Q->size1;
int i;
if (R->size1!=R->size2)return GSL_ENOTSQR;
if (R->size1!=n)return GSL_EINVAL;
if (Q->size1!=Q->size2)return GSL_ENOTSQR;
if (Q->size1!=n)return GSL_EINVAL;
gsl_vector* xi=gsl_vector_alloc(n);
gsl_vector* ei=gsl_vector_alloc(n);

for(i=0;i<n;i++){
gsl_vector_set_zero(ei);
gsl_vector_set(ei,i,1);/* constructing i'th unit vector*/
int status=QRsolve(Q,R,ei,xi);
fprintf(stderr,"%i\n",status);
gsl_matrix_set_col(B,i,xi);
}
gsl_vector_free(ei);
gsl_vector_free(xi);
    return GSL_SUCCESS;
}