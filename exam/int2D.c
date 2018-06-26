#include"math.h"
#include"stdio.h"
#include"assert.h"
double unity(double x, double y,double func(double x,double y)){return func(x,y);}
double anotb (double x, double y,double func(double x,double y)){return func(x/(1+x),y)/pow(1+x,2);}
double bnota(double x, double y,double func(double x,double y)){return func(x/(1-x),y)/pow(1-x,2);}
double doubleinf(double x, double y,double func(double x,double y)){return func(x/(1-x*x),y)*(1+x*x)/pow(1-x*x,2);}

double adaptive_slave(double func(double x,double y),double transform(double x, double y, double func(double x,double y)),double y,
        double a,double b,double f2,double f3,double acc, double eps,int current){

assert(current<100000);

double x1=a+(b-a)/6;
double x4=b-(b-a)/6;
double f1=transform(x1,y,func),f4=transform(x4,y,func);
double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
double q=(f1+f2+f3+f4)/4*(b-a);
double error=fabs(Q-q);
double tolerance=acc+eps*fabs(Q);
if(error<tolerance) return Q;
else{
double Q1=adaptive_slave(func,transform,y,a,a+(b-a)/2,f1,f2,acc/sqrt(2),eps,current+1);
double Q2=adaptive_slave(func,transform,y,b-(b-a)/2,b,f3,f4,acc/sqrt(2),eps,current+1);
return Q1+Q2;
}
}

double integrator(double func(double x,double y), double y,double a,double b,double acc, double eps){
double Q,Q1,Q2,f1,f2;
if (b==a) return 0;
if (b<a){Q=integrator(func,y,b,a,acc,eps);
return -Q;}
else{
int ainf=isinf(a);
int binf=isinf(b);

if (ainf==0){if (binf==0){
f1=func(a+(b-a)/3,y);
f2=func(b-(b-a)/3,y);
Q=adaptive_slave(func,unity,y,a,b,f1,f2,acc,eps,0);}
else{ 
    Q1=integrator(func,y,a,0,acc,eps);
    f1=bnota(1.0/3,y,func);
    f2=bnota(2.0/3,y,func);
    Q2=adaptive_slave(func,bnota,y,0,1,f1,f2,acc,eps,0);
    Q=Q1+Q2;}
}
else{if (binf==0) {Q1=integrator(func,y,0,b,acc,eps);
    f1=anotb(-1.0/3,y,func);
    f2=anotb(-2.0/3,y,func);
    Q2=adaptive_slave(func,anotb,y,-1,0,f1,f2,acc,eps,0);
    Q=Q1+Q2;}
    
   else {f1=doubleinf(-1.0/3,y,func); f2=doubleinf(1.0/3,y,func);
       Q=adaptive_slave(func,doubleinf,y,-1,1,f1,f2,acc,eps,0);
}}

  
    return Q;
}}

double slave2D(double func(double x,double y),double a,double b,double d(double y),double u(double y),
    double f2,double f3,double acc, double eps,int current){

assert(current<100000);
 
double y1=a+(b-a)/6;
double y4=b-(b-a)/6;
double f1=integrator(func,y1,d(y1),u(y1),acc,eps);
double f4=integrator(func,y4,d(y4),u(y4),acc,eps);
double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
double q=(f1+f2+f3+f4)/4*(b-a);
double error=fabs(Q-q);
double tolerance=acc/pow(sqrt(2),current)+eps*fabs(Q);/*Adjusting accuracy to take accuont for increased number of absiccas*/
if(error<tolerance) return Q;
else{
double Q1=slave2D(func,a,a+(b-a)/2,d,u,f1,f2,acc,eps,current+1);
double Q2=slave2D(func,b-(b-a)/2,b,d,u,f3,f4,acc,eps,current+1);
return Q1+Q2;
}
}

double int2D(double func(double x, double y), double a, double b, double d(double y), double u(double y),double acc,double eps){
if (b==a) return 0;
double Q;
if (b<a){Q=int2D(func,b,a,d,u,acc,eps);
return -Q;}
else{
assert(isinf(a)==0);
assert(isinf(b)==0);

double y2=a+(b-a)/3;
double y3=a+2*(b-a)/3;

double f2=integrator(func,y2,d(y2),u(y2),acc,eps);
double f3=integrator(func,y3,d(y3),u(y3),acc,eps);

Q=slave2D(func, a,b,d,u,f2,f3,acc,eps,0);
return Q;}
}
