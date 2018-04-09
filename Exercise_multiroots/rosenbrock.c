#include"stdlib.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_multiroots.h"
#include"stdio.h"
#include"gsl/gsl_errno.h"

int rosenbrock_gradient(gsl_vector* v, void* params, gsl_vector* f){
double x=gsl_vector_get(v,0), y=gsl_vector_get(v,1);
gsl_vector_set(f,0,-2*(1-x)-400*x*(y-x*x));
gsl_vector_set(f,1,200*(y-x*x));
return GSL_SUCCESS;
}

void findmin(gsl_vector* result){
int flag;
int iter=0,i_max=1000;
const gsl_multiroot_fsolver_type* T=gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver* S =gsl_multiroot_fsolver_alloc(T,2);
gsl_multiroot_function F;

F.f=(&rosenbrock_gradient);
F.n=2;
F.params=NULL;

gsl_vector* guess=gsl_vector_alloc(2);
gsl_vector_set(guess,0,2.0);
gsl_vector_set(guess,1,-2.0);

gsl_multiroot_fsolver_set(S,&F,guess);
do{iter ++;
flag=gsl_multiroot_fsolver_iterate(S);
if (flag){printf("Iteration error\n");break;}
flag=gsl_multiroot_test_residual(S->f,1e-9);
if (flag==GSL_SUCCESS){break;}
printf("%g %g %g %g\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1),gsl_vector_get(S->f,0),gsl_vector_get(S->f,1));
}while(flag==GSL_CONTINUE || iter<i_max);
if (iter==i_max) printf("Iteration limit reached\n");
/*
fprintf(stderr,"Status =%s",gsl_strerror(flag));*/

double x0=gsl_vector_get(S->x,0);
double y0=gsl_vector_get(S->x,1);
gsl_vector_set(result,0,x0);
gsl_vector_set(result,1,y0);

gsl_multiroot_fsolver_free(S);
gsl_vector_free(guess);
printf("\n\n");
printf("Number of iterations: %i\n",iter);
}
