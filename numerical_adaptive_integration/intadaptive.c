#include"math.h"
#include"stdio.h"
#include"assert.h"
double unity(double x, double func(double)){return func(x);}


double adaptive_slave(double func(double),double transform(double),double a,double b,double f2,double f3,double acc, double eps,void* iter,int current){

assert(current<100000);
 *((int*)iter)=current+1;
double x1=a+(b-a)/6;
double x4=b-(b-a)/6;
double f1=transform(x,func),f4=transform(x,func);
double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
double q=(f1+f2+f3+f4)/4*(b-a);
double error=fabs(Q-q);
double tolerance=acc+eps*fabs(Q);
if(error<tolerance) return Q;
else{
double Q1=adaptive_slave(func,transform,a,a+(b-a)/2,f1,f2,acc/sqrt(2),eps,iter, current+1);
double Q2=adaptive_slave(func,transform,b-(b-a)/2,b,f3,f4,acc/sqrt(2),eps,iter, current+1);
return Q1+Q2;
}
}

double integrator(double func(double),double a,double b,double acc, double eps, int printboolean){
int ainf=isinf(a);
int binf=isinf(b);
int iter=0;
assert(a<b);
if (ainf==0){if (binf==0)double(*trans)(double)=&unity;
else{ }

}
double f1=trans(a+(b-a)/3,func);
double f2=trans(b-(b-a)/3,func);
double Q=adaptive_slave(func,trans,a,b,f1,f2,acc,eps,(void*)&iter,0);
if(printboolean==1){ printf("Number of recursion levels in integration: %i.\n", iter);}
  
    return Q;
}
