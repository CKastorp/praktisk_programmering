#include"math.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_odeiv2.h"
#include"stdio.h"
int orbit_ode(double phi, double u[], double dudphi[],void* params){
double epsilon =*(double*)params;
dudphi[0]=u[1];
dudphi[1]=1-u[0]+epsilon*u[0]*u[0];
return GSL_SUCCESS;
}

double orbit_calculator(double x,double uprime, double epsilon){
double y[2];
y[0]=1;
y[1]=uprime;
double t=0,acc=1e-6,eps=1e-6,hstart=0.01;
gsl_odeiv2_system sys;
sys.function=orbit_ode;
sys.jacobian=NULL;
sys.dimension=2;
sys.params=(void*)&epsilon;

gsl_odeiv2_driver* driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,acc,eps);
gsl_odeiv2_driver_apply(driver,&t,x,y);
gsl_odeiv2_driver_free(driver);
return y[0];
}

int orbit(){
double x,y,r,phi;
for(phi=0;phi<2*M_PI;phi+=0.01){
r=1.0/orbit_calculator(phi,0,0);
x=r*cos(phi);y=r*sin(phi);
printf("%g %g\n",x,y);
}

printf("\n\n");
for(phi=0;phi<2*M_PI;phi+=0.01){
r=1.0/orbit_calculator(phi,-0.5,0);
x=r*cos(phi);y=r*sin(phi);
printf("%g %g\n",x,y);
}
printf("\n\n");

for (phi=0;phi<30;phi+=0.01){
r=1.0/orbit_calculator(phi,-0.5,0.01);x=r*cos(phi);
y=r*sin(phi);
printf("%g %g\n",x,y);
}

return GSL_SUCCESS;
}
