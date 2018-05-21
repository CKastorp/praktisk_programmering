#include"gsl/gsl_errno.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_matrix.h"
#include"math.h"
#include"stdio.h"

int newton(void func(gsl_vector* x, gsl_vector* fx),gsl_vector* x0, double dx,double epsilon);
void testfunction(gsl_vector* v,gsl_vector*fx){
double x=gsl_vector_get(v,0);
double y=gsl_vector_get(v,1);
double A=10000.0;
gsl_vector_set(fx,0,A*x*y-1);
gsl_vector_set(fx,1,exp(-x)+exp(-y)-1-1.0/A);
}

void rosenbrock_gradient(gsl_vector*v,gsl_vector* fx){
double x=gsl_vector_get(v,0),y=gsl_vector_get(v,1);
gsl_vector_set(fx,0,2-2*x-400*x*(y-x*x));
gsl_vector_set(fx,1,200*(y-x*x));
}

void himmelblau_gradient(gsl_vector* v,gsl_vector* fx){
    double x=gsl_vector_get(v,0),y=gsl_vector_get(v,1);
    double dx=4*x*(x*x+y-11)+2*(x+y*y-7);
    double dy=2*(x*x+y-11)+4*y*(x+y*y-7);
    gsl_vector_set(fx,0,dx); gsl_vector_set(fx,1,dy);
}

int main(){
int d=2;
gsl_vector* x=gsl_vector_alloc(d);
gsl_vector* fx=gsl_vector_alloc(d);

gsl_vector_set(x,0,2);
gsl_vector_set(x,1,-2);
double dx=0.000001;
double epsilon=0.00001;

int status=newton(testfunction,x,dx,epsilon);
testfunction(x,fx);
if(status==11) fprintf(stderr,"Newton root finding not converged, equations.\n");
printf("Solution to equations: f(x)=(%g,%g) at x=(%g,%g).\n",gsl_vector_get(fx,0),gsl_vector_get(fx,1),gsl_vector_get(x,0),gsl_vector_get(x,1));

gsl_vector_set(x,0,1.5);
gsl_vector_set(x,1,0.2);
status=newton(rosenbrock_gradient,x,dx,epsilon);
if(status==11) fprintf(stderr,"Newton root finding not converged, Rosenbrock.\n");
printf("Rosenbrock minimum: (%g,%g). Analytic result: (1,1).\n",gsl_vector_get(x,0),gsl_vector_get(x,1));

gsl_vector_set(x,0,-3);
gsl_vector_set(x,1,-3);
status=newton(himmelblau_gradient,x,dx,epsilon);
if(status==11) fprintf(stderr,"Newton root finding not converged, Himmelblau.\n");
printf("Himmelblau minimum found: (%g,%g).\n",gsl_vector_get(x,0),gsl_vector_get(x,1));

gsl_vector_free(x);
gsl_vector_free(fx);
    return 0;
}