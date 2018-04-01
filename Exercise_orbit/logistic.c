#include"stdio.h"
#include"math.h"
#include"gsl/gsl_odeiv2.h"
#include"gsl/gsl_errno.h"
int logistic_ode(double x, double y[],double dydx[],void* params){
dydx[0]=y[0]*(1-y[0]);
return GSL_SUCCESS;
}

double logistic_calculator(double x){
double y[2]={0.5,0.25};
double t=0,acc=1e-6,eps=1e-6,hstart=0.01;
gsl_odeiv2_system sys;
sys.function=logistic_ode;
sys.jacobian=NULL;
sys.dimension=1;
sys.params=NULL;

gsl_odeiv2_driver* driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,acc,eps);
gsl_odeiv2_driver_apply(driver,&t,x,y);
gsl_odeiv2_driver_free(driver);
return y[0];
}

int logistic(){
double x;
for(x=0;x<=3;x+=0.01){
printf("%g %g\n",x,logistic_calculator(x));
}
return GSL_SUCCESS;
}
