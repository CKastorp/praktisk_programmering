#include"stdio.h"
#include"math.h"
#include"stdlib.h"
int binsearch(double x[], double z, int size);

double lininterp_integrator(double x[], double y[], double z,int n){

double result, deltaX, deltaY, A, y0=0;
int index, i;

index=binsearch(x,z,n);
if(index>0){
for(i=0;i<index;i++){
deltaX=x[i+1]-x[i];
deltaY=y[i+1]-y[i];
y0+=deltaX*deltaY/2+deltaX*y[i];
}}

deltaX=x[index+1]-x[index];
deltaY=y[index+1]-y[index];

A=deltaY/deltaX;


result=y0+A*(z-x[index])*(z-x[index])/2+y[index]*(z-x[index]);

return result;
}
