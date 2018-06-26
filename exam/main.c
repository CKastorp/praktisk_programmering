#include"math.h"
#include"stdio.h"
double testfunc(double x, double y){return 1.0;}
double dummyd(double x){return 0;}
double dummyu(double x){return 1;}
double lowercircle(double x){return -sqrt(1-x*x);}
double uppercircle(double x){return sqrt(1-x*x);}
double gauss2D(double x,double y){return exp(-x*x-y*y);}

double integrator(double func(double),double a,double b,double abs, double eps, int print);
double int2D(double func(double x, double y), double a, double b, double d(double x), double u(double x),double acc,double eps);

int main(){
    
    double eps=0.000001, abs=0.000001, a=0.0, b=1.0;
    double result=int2D(testfunc,a,b,dummyd,dummyu,abs,eps);
    printf("First test: I=%g. Should be 1.\n",result);
    a=-1.0;
    result=int2D(testfunc,a,b,lowercircle,uppercircle,abs,eps);
    printf("Area of unit disc: I=%g. Should be pi.\n",result);
    result=int2D(gauss2D,a,b,lowercircle,uppercircle,abs,eps);
    printf("2D Gauss integrated over a unit disc: I=%g.\n",result);

    return 0;}