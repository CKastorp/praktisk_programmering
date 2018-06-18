#include"math.h"
#include"stdio.h"
#include"assert.h"

double adaptive_slave(double func(double),double a,double b,double f2,double f3,double acc, double eps,void* iter,int current){


assert(current<100000);
 *((int*)iter)=current+1;
double x1=a+(b-a)/6;
double x4=b-(b-a)/6;
double f1=func(x1),f4=func(x4);
double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
double q=(f1+f2+f3+f4)/4*(b-a);
double error=fabs(Q-q);
double tolerance=acc+eps*fabs(Q);
if(error<tolerance) return Q;
else{
double Q1=adaptive_slave(func,a,a+(b-a)/2,f1,f2,acc/sqrt(2),eps,iter, current+1);
double Q2=adaptive_slave(func,b-(b-a)/2,b,f3,f4,acc/sqrt(2),eps,iter, current+1);
return Q1+Q2;
}
}

double integrator(double func(double),double a,double b,double acc, double eps, int printboolean){

int iter=0;

double f1=func(a+(b-a)/3);
double f2=func(b-(b-a)/3);
double Q=adaptive_slave(func,a,b,f1,f2,acc,eps,(void*)&iter,0);
if(printboolean==1){ printf("Number of recursion levels in integration: %i.\n", iter);}
  
    return Q;
}