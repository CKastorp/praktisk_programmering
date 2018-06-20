#include"math.h"
#include"stdio.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"

int sgn(double x){if (x<0)return -1; if (x==0) return 0; if (x>0)return 1;}
void rk45_step(double t,double h,gsl_vector* yt,gsl_vector* yth, void func(double t, gsl_vector* y,gsl_vector* dydt),gsl_vector* err){
int i,n=yt->size;
gsl_vector* temp=gsl_vector_alloc(n);
gsl_vector* k0=gsl_vector_alloc(n);
gsl_vector* k1=gsl_vector_alloc(n);
gsl_vector* k2=gsl_vector_alloc(n);
gsl_vector* k3=gsl_vector_alloc(n);
gsl_vector* k4=gsl_vector_alloc(n);
gsl_vector* k5=gsl_vector_alloc(n);

func(t,yt,k0);
for(i=0;i<n;i++)gsl_vector_set(temp,i,gsl_vector_get(yt,i)+0.5*gsl_vector_get(k0,i)*h);

func(t+0.25*h,temp,k1);
for(i=0;i<n;i++)gsl_vector_set(temp,i, gsl_vector_get(yt,i)+(3.0/32*gsl_vector_get(k0,i)+9.0/32*gsl_vector_get(k1,i))*h);

func(t+0.375*h,temp,k2);
for(i=0;i<<n;i++)gsl_vector_set(temp,i,gsl_vector_get(yt,i)+h*(1932.0/2197*gsl_vector_get(k0,i)
    -7200/2197*gsl_vector_get(k1,i)+7296/2197*gsl_vector_get(k2,i)));

func(t+12.0/13*h,temp,k3);
for(i=0;i<n;i++)gsl_vector_set(temp,i,gsl_vector_get(yt,i)+h*(439.0/216*gsl_vector_get(k0,i)-8*gsl_vector_get(k1,i)+
3680.0/513*gsl_vector_get(k2,i)-845.0/4104*gsl_vector_get(k3,i)));

func(t+h,temp,k4);
for(i=0;i<n;i++)gsl_vector_set(temp,i,gsl_vector_get(yt,i)+h*(-8.0/27*gsl_vector_get(k0,i)+2*gsl_vector_get(k1,i)-
3544.0/2565*gsl_vector_get(k2,i)+1859.0/4104*gsl_vector_get(k3,i)-11.0/40*gsl_vector_get(k4,i)));

func(t+0.5*h,temp,k5);
for(i=0;i<n;i++){gsl_vector_set(yth,i,gsl_vector_get(yt,i)+h*(16.0/135*gsl_vector_get(k0,i)+6656.0/12825*gsl_vector_get(k2,i)+
    28561.0/56430*gsl_vector_get(k3,i)-9.0/50*gsl_vector_get(k4,i)+2.0/55*gsl_vector_get(k5,i)));
gsl_vector_set(temp,i,gsl_vector_get(yt,i)+h*(25.0/216*gsl_vector_get(k0,i)+1408.0/2565*gsl_vector_get(k2,i)
    +2197.0/4104*gsl_vector_get(k3,i)-0.2*gsl_vector_get(k4,i)));
}

gsl_vector_memcpy(err,yth);
gsl_vector_sub(err,temp);

gsl_vector_free(k0);
gsl_vector_free(k1);
gsl_vector_free(k2);
gsl_vector_free(k3);
gsl_vector_free(k4);
gsl_vector_free(k5);
gsl_vector_free(temp);}

int rk45_driver(double t,double b,double h, gsl_vector* yt,void func(double t, gsl_vector* y,gsl_vector* dydt),double acc,double eps){
if (sgn(h)!=sgn(b-t)) h*=-1;
if(t>=b)return GSL_EINVAL;
int n=yt->size;
gsl_vector* err=gsl_vector_alloc(n);
gsl_vector* yth=gsl_vector_alloc(n);
int iter=0,quit=0,limit=100000000;
double ynorm,enorm,stepsize;

while(iter<limit){
iter++;
if (t>=b) break;
stepsize=4*h;
if(t+h>b) stepsize=b-t;

ynorm=gsl_blas_dnrm2(yt);
enorm=acc+ynorm*eps*2;/*Make sure that the error is large enough for the loop to initiate.*/

while(sgn(acc+ynorm*eps-enorm)<0){
rk45_step(t,stepsize,yt,yth,func,err);
if(stepsize<0.00000001){quit=1;break;}
stepsize*=0.5;
enorm=gsl_blas_dnrm2(err);

}
stepsize*=2;
if (quit==1)break;
h=stepsize;
gsl_vector_memcpy(yt,yth);
t+=stepsize;
/*fprintf(stderr,"t=%g, Y=(%g,%g). Step: %g.\n",t,gsl_vector_get(yt,0),gsl_vector_get(yt,1),stepsize);*/

}


gsl_vector_free(err);
gsl_vector_free(yth);
if (quit==1)return GSL_ETOL;
else{if(iter==limit) return GSL_EMAXITER;
  else  return GSL_SUCCESS;
}}