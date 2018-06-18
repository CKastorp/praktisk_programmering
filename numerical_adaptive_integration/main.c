#include"math.h"
#include"stdio.h"
double testfunc(double x){return 1.0;}
double squareroot(double x){return sqrt(x);}
double inversesquareroot(double x){return 1.0/sqrt(x);}
double logroot(double x){return log(x)/sqrt(x);}
double calcpi(double x){return 4*sqrt(1-(1-x)*(1-x));}
double integrator(double func(double),double a,double b,double abs, double eps, int print);

int main(){
    int print=1;
    double eps=0.000001, abs=0.000001, a=0.0, b=1.0;
    double result=integrator(testfunc,a,b,abs,eps,print);
    printf("First test: I=%g. Should be 1.\n",result);
    result=integrator(squareroot,a,b,abs,eps,1);
    printf("int_0^1 sqrt(x) dx= %g. Should be 2/3.\n",result);
    
    result=integrator(inversesquareroot,a,b,abs,eps,1);
    printf("int_0^1 1/sqrt(x) dx= %g. Should be 2.\n",result);
    result=integrator(logroot,a,b,abs,eps,1);
    printf("int_0^1 ln(x)/sqrt(x) dx= %g. Should be -4.\n",result);
    abs=0.000000001;eps=0.000000001;
    result=integrator(calcpi,a,b,abs,eps,1);
    printf("pi=%g.\n",result);
    return 0;}