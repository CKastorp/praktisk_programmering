#include"gsl/gsl_odeiv2.h"
#include"gsl/gsl_errno.h"
#include"math.h"
#include"stdio.h"
#include"assert.h"
#include"gsl/gsl_multiroots.h"

int ode_Swave(double r, double y[],double dydr[],void* params){
double e=*(double*)params;
dydr[0]=y[1];
dydr[1]=2*(-e-1/r)*y[0];

return GSL_SUCCESS;
}

double hydrogen(double r, double e){
assert(r>=0);
double t=1e-3;
if (r<t) return r-r*r;
double startstep=1e-3, abs=1e-6,rel=1e-6;

gsl_odeiv2_system sys;
sys.function =ode_Swave;
sys.jacobian=NULL;
sys.dimension=2;
sys.params=(void*)&e;
gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, startstep,abs,rel);
double y[]={t-t*t,1-2*t};

int flag=gsl_odeiv2_driver_apply(driver, &t, r,y);
if (flag!= GSL_SUCCESS) fprintf(stderr,"Error, odeiv:stutus=%s\n",gsl_strerror(flag));
gsl_odeiv2_driver_free(driver);
return y[0];
}

int hydrogen_set_function(gsl_vector* x,void* params, gsl_vector* f){
double r=*(double*)params;
double e=gsl_vector_get(x,0);
gsl_vector_set(f,0,hydrogen(r,e));
return GSL_SUCCESS;
}
int hydrogen_master(double rmax){

gsl_multiroot_fsolver* solver=gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,1);
gsl_multiroot_function F;
F.f=hydrogen_set_function;
F.n=1;
F.params=(void*)&rmax;

gsl_vector* x=gsl_vector_alloc(1);
gsl_vector_set(x,0,-0.8);

gsl_multiroot_fsolver_set(solver, &F, x);

int iter=0;
int flag;

do {
iter ++;
flag=gsl_multiroot_fsolver_iterate(solver);
if (flag!=GSL_SUCCESS){break;}
flag=gsl_multiroot_test_residual(solver->f,1e-9);
if (flag==GSL_SUCCESS){break;}
}while(iter<1000 && flag== GSL_CONTINUE);

double e=gsl_vector_get((solver->x),0);
printf("rmax: %g, E=%g\n",rmax,e);
printf("\n\n");

double r;
for(r=0;r<=rmax;r+=0.01)printf("%g %g %g\n",r,hydrogen(r,e),r*exp(-r));

gsl_vector_free(x);
gsl_multiroot_fsolver_free(solver);
return GSL_SUCCESS;
}
