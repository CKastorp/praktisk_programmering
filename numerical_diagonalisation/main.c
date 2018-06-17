#include"stdio.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"
void print_vector(gsl_vector* v);
void make_sym_matrix(gsl_matrix*A);
void print_matrix(gsl_matrix* A);
int eig_sweep(gsl_matrix* A, gsl_vector* values, gsl_matrix* V);

int main(int argc, char** argv){
int n;
if (argc>1) n=atof((argv[1]));
else return GSL_EINVAL;

gsl_matrix* A=gsl_matrix_alloc(n,n);
gsl_matrix* Acpy=gsl_matrix_alloc(n,n);
gsl_matrix* V=gsl_matrix_alloc(n,n);
gsl_vector* values=gsl_vector_alloc(n);
gsl_matrix* VTAV=gsl_matrix_alloc(n,n);
make_sym_matrix(A);
printf("Matrix A:\n");
print_matrix(A);
gsl_matrix_memcpy(Acpy,A);

int status=eig_sweep(A,values,V);

printf("Eigenvalues after procedure:\n");print_vector(values);
printf("Matrix V after procedure:\n");print_matrix(V);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,Acpy,0,A);
/*printf("Matrix product V^TA:\n");print_matrix(VTA);*/
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,V,0,VTAV);
printf("Matrix product V^TAV:\n");print_matrix(VTAV);

gsl_matrix_free(VTAV);
gsl_vector_free(values);
gsl_matrix_free(V);
gsl_matrix_free(Acpy);
gsl_matrix_free(A);
    return 0;
}