#include"gsl/gsl_odeiv2.h"
#include"gsl/gsl_errno.h"
#include"math.h"
#include"stdio.h"
int exponential_ode(double x,double y[],double dydx[],void* params){
dydx[0]=y[0];
dydx[1]=y[0];
return GSL_SUCCESS;
}

double exponential_calculator(double x){
if(x==0) return 1;
if(x<0) return 1.0/exponential_calculator(-x);
if (x>1) {double Q=exponential_calculator(x/2); return Q*Q;}
double y[2]={1,1};
double t=0,acc=1e-6,eps=1e-6,hstart=0.01;
gsl_odeiv2_system sys;
sys.function=exponential_ode;
sys.jacobian=NULL;
sys.dimension=2;
sys.params=NULL;

gsl_odeiv2_driver* driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,acc,eps);
gsl_odeiv2_driver_apply(driver,&t,x,y);
gsl_odeiv2_driver_free(driver);
return y[0];
}

int main(){
double x;
for(x=-2;x<=3;x+=0.2){
printf("%g %g %g\n",x,exponential_calculator(x),exp(x));
}
return 0;
}
