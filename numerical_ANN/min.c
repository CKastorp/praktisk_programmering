#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"math.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"

double dotproduct(gsl_vector* A,gsl_vector* B){double result=0;
    for(int i=0;i<A->size;i++){
    double Ai=gsl_vector_get(A,i);
    double Bi=gsl_vector_get(B,i);
result+=Ai*Bi;
}
return result;
}

void num_grad(gsl_vector* x,gsl_vector* g,double f(gsl_vector* x)){
double dx=0.00000015;
double fx=f(x);
for (int i=0;i<x->size;i++){
double xi=gsl_vector_get(x,i);
gsl_vector_set(x,i,xi+dx);
gsl_vector_set(g,i,(f(x)-fx)/dx);
gsl_vector_set(x,i,xi);
}
/*fprintf(stderr,"Gradient:(%g,%g).\n",gsl_vector_get(g,0),gsl_vector_get(g,1));*/
}

double linsearch(double func(gsl_vector* x),gsl_vector*x,gsl_vector*s,gsl_vector*df,double lambda_min){
double lambda=1.0;
double alpha=0.01;
/*int n=x->size;*/
while(lambda>lambda_min){
double RHS=func(x)+alpha*dotproduct(s,df);
gsl_vector_add(s,x);
double LHS=func(s);
gsl_vector_sub(s,x);
if (LHS < RHS)break;
gsl_vector_scale(s,0.5);
lambda*=0.5;
}
return lambda;
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

 /*   printf("Quasi-newtonian, number of steps: %i.\n", steps);*/
if(steps<max_step) return GSL_SUCCESS;
    else return GSL_EMAXITER;
}
