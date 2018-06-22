#include"gsl/gsl_blas.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_errno.h"
#include"math.h"

void print_matrix(gsl_matrix* M);
void print_vector(gsl_vector* x);
int QRsolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
int QRfactor(gsl_matrix* A,gsl_matrix* R);

int newton(void func(gsl_vector* x, gsl_vector* fx),gsl_vector* x, double dx,double epsilon){
int d=x->size;
gsl_vector* fx=gsl_vector_alloc(d);
gsl_matrix* J=gsl_matrix_alloc(d,d);
gsl_matrix* R=gsl_matrix_alloc(d,d);
gsl_vector* df=gsl_vector_alloc(d);
gsl_vector* z=gsl_vector_alloc(d);
gsl_vector* fz=gsl_vector_alloc(d);
gsl_vector* step=gsl_vector_alloc(d);
gsl_vector* Dx=gsl_vector_alloc(d);

double lambda;
int iter=0, iter_max=10000,i,j,n=0;
while(iter<iter_max){
iter++;
func(x,fx);

for (j=0;j<d;j++){
    double xj=gsl_vector_get(x,j);
gsl_vector_set(x,j,xj+dx);
func(x,df);
gsl_vector_sub(df,fx);
for(i=0;i<d;i++)gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx);
gsl_vector_set(x,j,xj);
}

QRfactor(J,R);

gsl_vector_scale(fx,-1);
int status=QRsolve(J,R,fx,Dx);

gsl_vector_scale(fx,-1);
double fx_abs=gsl_blas_dnrm2(fx);
n=0;
while(n<7){
    lambda=1.0/pow(2,n);
gsl_vector_memcpy(z,x);
gsl_vector_memcpy(step,Dx);
gsl_vector_scale(step,lambda);
gsl_vector_add(z,step);
func(z,fz);
if(gsl_blas_dnrm2(fz)<fx_abs*(1-lambda/2)) break;
n++;
}

gsl_vector_add(x,step);
gsl_vector_memcpy(fx,fz);



if(gsl_blas_dnrm2(step)<dx||gsl_blas_dnrm2(fx)<epsilon) break;
}
printf("Iterations: %i.\n",iter);
gsl_vector_free(Dx);
gsl_vector_free(step);
gsl_vector_free(fz);
gsl_vector_free(z);
gsl_vector_free(df);
gsl_matrix_free(R);
gsl_matrix_free(J);
gsl_vector_free(fx);

if (iter<iter_max) return GSL_SUCCESS;
else return GSL_EMAXITER;
}

