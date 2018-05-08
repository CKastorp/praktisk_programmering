#include"stdio.h"
#include"math.h"
int binsearch(double x[], double z,int n);

double linear_interpolation(double x[], double y[], double z, int n){

double result, deltaX, deltaY, A;
int index;

index=binsearch(x,z,n);
deltaX=x[index+1]-x[index];
deltaY=y[index+1]-y[index];
A=deltaY/deltaX;

result=A*(z-x[index])+y[index];

return result;
}
