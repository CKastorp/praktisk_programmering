#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"math.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"
void print_matrix(gsl_matrix*A);
double dotproduct(gsl_vector* A,gsl_vector* B);
void num_grad(gsl_vector* x,gsl_vector* g,double f(gsl_vector* x));
double linsearch(double func(gsl_vector* x),gsl_vector*x,gsl_vector*s,gsl_vector*df,double lambda_min);
int QRsolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
int QRfactor(gsl_matrix* A,gsl_matrix* R);

int newton(double func(gsl_vector* x), void set_hessian(gsl_vector*x,gsl_matrix*H),gsl_vector* x0,gsl_vector* df,double epsilon){
int n=x0->size;
gsl_matrix* H=gsl_matrix_alloc(n,n);
gsl_matrix* R=gsl_matrix_alloc(n,n);
gsl_vector* s=gsl_vector_alloc(n);
int steps=0,limit=10000;
double lambda_min=1.0/1000000;
num_grad(x0,df,func);

while(steps<limit){
steps++;
set_hessian(x0,H);

QRfactor(H,R);
QRsolve(H,R,df,s);
gsl_vector_scale(s,-1);
double lambda=linsearch(func,x0,s,df,lambda_min);
gsl_vector_add(x0,s);
num_grad(x0,df,func);
if(sqrt(dotproduct(df,df))<epsilon)break;
}

printf("Newtonian optimisation, number of steps: %i.\n",steps);
gsl_matrix_free(R);
gsl_vector_free(s);
gsl_matrix_free(H);
if(steps<limit) return GSL_SUCCESS;
    else return GSL_EMAXITER;
}


int qnewton(double func(gsl_vector*x),gsl_vector* x0,gsl_vector* df,double epsilon){
int n=x0->size;
int steps=0, max_step=100000;
double lambda_min=1.0/1000000000;
gsl_matrix* B=gsl_matrix_alloc(n,n);
gsl_vector* s=gsl_vector_alloc(n);
gsl_vector* y=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);

gsl_matrix_set_identity(B);

while(steps<max_step){
    steps++;
num_grad(x0,df,func);
  /*  fprintf(stderr,"Position: (%g,%g). Gradient: (%g,%g).\n",gsl_vector_get(x0,0),gsl_vector_get(x0,1),gsl_vector_get(df,0),gsl_vector_get(df,1));*/
gsl_blas_dgemv(CblasNoTrans,-1,B,df,0,s);
double lambda=linsearch(func,x0,s,df,lambda_min);
if (lambda<=lambda_min)gsl_matrix_set_identity(B);
gsl_vector_add(x0,s);
num_grad(x0,y,func);

    if(gsl_blas_dnrm2(y)<epsilon)break;
gsl_vector_sub(y,df);/*y=grad(x+dx)-grad(x)*/

gsl_vector_memcpy(u,s);
gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);

double product=dotproduct(u,y);
if(fabs(product)>0.00000000001){


gsl_blas_dger(1.0/product,u,u,B);

}
}
gsl_vector_free(u);
gsl_vector_free(y);
gsl_vector_free(s);
gsl_matrix_free(B);

    printf("Quasi-newtonian, number of steps: %i.\n", steps);
if(steps<max_step) return GSL_SUCCESS;
    else return GSL_EMAXITER;
}
