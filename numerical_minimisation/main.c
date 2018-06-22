#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"math.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"
void rosen_hessian(gsl_vector*x0,gsl_matrix* H);
double rosen(gsl_vector* x);
void himmel_hessian(gsl_vector*v,gsl_matrix*H);
double himmelblau(gsl_vector* v);
int newton(double func(gsl_vector* x), void set_hessian(gsl_vector*x,gsl_matrix*H),gsl_vector* x0,gsl_vector* df,double epsilon);/* For analytic Hessian matrices.*/
int qnewton(double func(gsl_vector*x),gsl_vector* x0,gsl_vector* df,double epsilon);

int main(){
    int n=2;
gsl_vector* xstart=gsl_vector_alloc(n);
gsl_vector* df=gsl_vector_alloc(n);
double epsilon=0.000001;
gsl_vector_set(xstart,0,-2);gsl_vector_set(xstart,1,2);

int status=newton(rosen,rosen_hessian,xstart,df,epsilon);
if (status==11)fprintf(stderr,"Rosenbrock minimum not converged.\n");
printf("Solution found to Rosenbrock minimum: (%g,%g). Analytic: (1,1).\n",gsl_vector_get(xstart,0),gsl_vector_get(xstart,1));

gsl_vector_set(xstart,0,5); gsl_vector_set(xstart,1,5);

status=newton(himmelblau,himmel_hessian,xstart,df,epsilon);
if (status==11)fprintf(stderr,"Himmelblau minimum not converged.\n");
printf("Solution found to Himmelblau minimum: (%g,%g).\n",gsl_vector_get(xstart,0),gsl_vector_get(xstart,1));

gsl_vector_set(xstart,0,1.5);gsl_vector_set(xstart,1,0.2);
epsilon=0.000001;
status=qnewton(rosen,xstart,df,epsilon);
if (status==11)fprintf(stderr,"Rosenbrock minimum not converged (qnewton).\n");
printf("Solution found to Rosenbrock minimum (Quasi-Newton): (%g,%g). Analytic: (1,1).\n",gsl_vector_get(xstart,0),gsl_vector_get(xstart,1));

gsl_vector_set(xstart,0,5);gsl_vector_set(xstart,1,5);
status=qnewton(himmelblau,xstart,df,epsilon);
if (status==11)fprintf(stderr,"Rosenbrock minimum not converged (qnewton).\n");
printf("Solution found to Himmelblau minimum (Quasi-Newton): (%g,%g).\n",gsl_vector_get(xstart,0),gsl_vector_get(xstart,1));
gsl_vector_free(df);
gsl_vector_free(xstart);
    return 0;
}