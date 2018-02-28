#include "stdio.h"
#include "limits.h"
#include "float.h"
#include "math.h"
void summation(void){
int n=1;
int Nmax=INT_MAX/2;
float sum_up_float=0;
for (n=1;n<=Nmax;n++){sum_up_float+=1.0/n;}
float sum_down_float=0;
for (n=Nmax;n>=1;n--){sum_down_float+=1.0/n;}
printf("Float sum for increasing n: %f.\n Sum for decreasing n: %f.\n",sum_up_float,sum_down_float);
printf("Difference is due to rounding errors in float as 1/n approaches 0.\n");
printf("Due to finite sensitivity, the sum should converge.\n");
Nmax=2*(INT_MAX/3);
n=1;
sum_up_float=0;
while (n<=Nmax){sum_up_float+=1.0/n; n++;}
printf("Summation for increasing n, larger Nmax: %f.\n",sum_up_float);
Nmax=INT_MAX/2;
printf("Same procedure four doubles:\n");
double sum_up_double=0;
double sum_down_double=0;
for(n=1;n<=Nmax;n++){sum_up_double+=1.0/n;}
for(n=Nmax;n>=1;n--){sum_down_double+=1.0/n;}
printf("Increasing n: %g. Decreasing n: %g\n",sum_up_double,sum_down_double);
printf("The results are identical due to better precision in the 'double' type.\n");
}
