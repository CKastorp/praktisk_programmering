#include "stdio.h"
#include "math.h"
#include "complex.h"
int main() {

double complex q=csqrt(-2);
double complex z=cpow(M_E,I);
double complex v=cpow(M_E,I*M_PI);
double complex w=cpow(I,M_E);
printf("Gamma(5)=%g, j1(0.5)=%g\n",tgamma(5), j1(0.5));

printf("Funny complex numbers: %g%+gi, %g%+gi, %g%+gi, and %g%+gi\n",
creal(q),cimag(q),creal(z),cimag(z),creal(v),cimag(v),creal(w),cimag(w));

long double x=0.11111111111111111111111111111111L;
double y=x;
float yy=x;
printf("long %.25Lg, double %.25lg, float %.25g\n",x,y,yy);

return 0;
}
