#include"gsl/gsl_blas.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"math.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_matrix.h"
int QRsolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
int QRfactor(gsl_matrix* A,gsl_matrix* R);
void print_matrix(gsl_matrix*A);
void random_matrix(gsl_matrix* A);

int main(){
int n=4;
int m=3;
gsl_matrix* A=gsl_matrix_alloc(n,m);
gsl_matrix* R=gsl_matrix_alloc(m,m);
gsl_matrix* Q=gsl_matrix_alloc(n,m);
gsl_matrix* QTQ=gsl_matrix_alloc(m,m);
gsl_matrix* QRtest=gsl_matrix_alloc(n,m);
gsl_vector* x=gsl_vector_alloc(m);
gsl_vector* b=gsl_vector_alloc(n);

random_matrix(A);
gsl_matrix_memcpy(Q,A);
printf("Matrix A:\n");
print_matrix(A);
fprintf(stderr,"Doing QR\n");
int status=QRfactor(Q,R);

printf("Matrix Q after QR procedure:\n");
print_matrix(Q);
printf("Matrix R after QR procedure:\n");
print_matrix(R);

gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Q,Q,0,QTQ);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Q,R,0,QRtest);
printf("Matrix product Q^T Q:\n");
print_matrix(QTQ);
printf("Matrix product QR:\n");
print_matrix(QRtest);

fprintf(stderr,"Solving\n");
printf("Test for constructed b, x=(1,2,3):\n");
gsl_vector_set_zero(b);
for(int i=0;i<n;i++){
    double temp=0;
for(int j=0;j<m;j++){
temp+=gsl_matrix_get(A,i,j)*(j+1);
}
gsl_vector_set(b,i,temp);
}
QRsolve(Q,R,b,x);
printf("Solution found by routine: (%g,%g,%g).\n",gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2));
gsl_vector_free(x);
gsl_vector_free(b);
gsl_matrix_free(QRtest);
gsl_matrix_free(QTQ);
gsl_matrix_free(Q);
gsl_matrix_free(R);
gsl_matrix_free(A);
    return 0;
}