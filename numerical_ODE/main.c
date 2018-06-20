#include"gsl/gsl_vector.h"
#include"math.h"
#include"stdio.h"
#include"gsl/gsl_errno.h"
int rk45_driver(double t,double b,double h, gsl_vector* yt,void func(double t, gsl_vector* y,gsl_vector* dydt),double acc,double eps);

void mysine(double t,gsl_vector* y, gsl_vector* dydt){
gsl_vector_set(dydt,0,gsl_vector_get(y,1));
gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

void myairy(double t,gsl_vector* y,gsl_vector*dydt){
    gsl_vector_set(dydt,0,gsl_vector_get(y,1));
    gsl_vector_set(dydt,1,t*gsl_vector_get(y,0));
}

int main(){
int n=2;
gsl_vector* y=gsl_vector_alloc(n);

double t=0;

double h=0.02;
double acc=0.000001;
double eps=0.000001;

for(double k=1.0/16;k<15;k+=1.0/16){
gsl_vector_set(y,0,0);gsl_vector_set(y,1,1);
int status=rk45_driver(t,k,h,y,mysine,acc,eps);
printf("%g %g\n",k,gsl_vector_get(y,0));
}
printf("\n\n");
t=-15;
for(double k=-15+1.0/16;k<5;k+=1.0/16){
gsl_vector_set(y,0,0.3);gsl_vector_set(y,1,0);
int status=rk45_driver(t,k,h,y,myairy,acc,eps);
printf("%g %g\n",k,gsl_vector_get(y,0));}

gsl_vector_free(y);

    return 0;
}