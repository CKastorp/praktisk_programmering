#include"stdio.h"
#include"math.h"
#include"gsl/gsl_sf_bessel.h"
double linear_interpolation(double x[],double y[], double z, int n);
double lininterp_integrator(double x[], double y[], double z, int n);
int binsearch(double x[],double z, int size);

int main(){
int n=11;
double x[n],y[n];

for(int i=0;i<n;i++){
x[i]=i;
y[i]=i*i;
printf("%g %g\n",x[i],y[i]);
}
printf("\n\n");
fprintf(stderr,"Interpolation\n");
for(double k=0;k<n-1;k+=1.0/32)printf("%g %g %g\n",k,k*k,linear_interpolation(x,y,k,n));

fprintf(stderr,"Integration\n");
printf("\n\n");
for(double k=0;k<n-1;k+=1.0/32)printf("%g %g %g\n",k,k*k*k/3,lininterp_integrator(x,y,k,n));

return 0;
}
