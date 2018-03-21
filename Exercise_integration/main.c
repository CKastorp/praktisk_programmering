#include"stdio.h"
#include"stdlib.h"
#include"math.h"
#include"gsl/gsl_integration.h"
#include"gsl/gsl_errno.h"
double integration(void);
double dirac(double a);
int main(){
double result = integration();
printf("Result of integration is %g.\n",result);
printf("\n\n");
for (double i=0.1;i<3;i+=0.05)printf("%g %g\n",i,dirac(i));

return 0;
}
