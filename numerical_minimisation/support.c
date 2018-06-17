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
double alpha=0.001;
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
void print_matrix(gsl_matrix* Q){
int i,j;
for (i=0;i<Q->size1;i++){
for (j=0;j<Q->size2;j++){
printf("%g ",gsl_matrix_get(Q,i,j));
}
printf("\n");
}
printf("\n");
}
double matrix_column_dotproduct(gsl_matrix* A,int i,gsl_matrix* B,int j){
double result=0;
for(int n=0;n<A->size1;n++)result+=gsl_matrix_get(A,n,i)*gsl_matrix_get(B,n,j);

    return result;
}

double rosen(gsl_vector*v){
double x=gsl_vector_get(v,0);
double y=gsl_vector_get(v,1);
return pow(1-x,2)+100*pow(y-x*x,2);
}

void rosen_hessian(gsl_vector*v,gsl_matrix*H){
double x=gsl_vector_get(v,0);
double y=gsl_vector_get(v,1);

gsl_matrix_set(H,0,0,2.0-400.0*(y-x*x)+800.0*x*x);
gsl_matrix_set(H,0,1,-400.0*x);
gsl_matrix_set(H,1,0,-400.0*x);
gsl_matrix_set(H,1,1,200.0);
}

double himmelblau(gsl_vector* v){
double x=gsl_vector_get(v,0);
double y=gsl_vector_get(v,1);
double result= pow(x*x+y-11,2)+pow(x+y*y-7,2);
return result;
}

void himmel_hessian(gsl_vector*v,gsl_matrix*H){
double x=gsl_vector_get(v,0);
double y=gsl_vector_get(v,1);
gsl_matrix_set(H,0,0,8*x*x+2+4*(x*x+y-11));
gsl_matrix_set(H,0,1,4*x+4*y);
gsl_matrix_set(H,1,0,4*x+4*y);
gsl_matrix_set(H,1,1,2+8*y*y+4*(x+y*y-7));
}