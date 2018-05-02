#include"stdio.h"
#include"math.h"
#include"gsl/gsl_multimin.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_errno.h"

double rosenbrock(gsl_vector* v,void* params){
double x=gsl_vector_get(v,0), y=gsl_vector_get(v,1);
printf("%g %g \n",x,y);
return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}



int findmin(void){
int d=2, iter=0, status, iteration_status;
double acc=0.01;

gsl_multimin_fminimizer* system=gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,d);
gsl_multimin_function F;
F.f=rosenbrock;
F.n=d;

gsl_vector* start=gsl_vector_alloc(d);
gsl_vector_set(start,0,-2);
gsl_vector_set(start,1,2);

gsl_vector* stepsize=gsl_vector_alloc(d);
gsl_vector_set(stepsize,0,0.5);
gsl_vector_set(stepsize,1,0.5);

gsl_multimin_fminimizer_set(system,&F,start,stepsize);

while(iter<1000){
iter++;
iteration_status=gsl_multimin_fminimizer_iterate(system);
if (iteration_status!=0){fprintf(stderr,"Unable to improve\n");break;}
status=gsl_multimin_test_size(system->size,acc);
if(status ==GSL_SUCCESS){fprintf(stderr,"Converged\n"); break;}
}

double xmin=gsl_vector_get(system->x,0), ymin=gsl_vector_get(system->x,1);
printf("\n\n");
printf("Rosenbrock minimum: (%g, %g). Number of iterations: %i.\n",xmin,ymin,iter);
printf("Analytic result: (1,1).\n");

gsl_vector_free(stepsize);
gsl_vector_free(start);
gsl_multimin_fminimizer_free(system);
return GSL_SUCCESS;
}
