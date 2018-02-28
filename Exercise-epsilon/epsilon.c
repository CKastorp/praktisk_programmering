#include "stdio.h"
#include "limits.h"
#include "math.h"
#include "float.h"
void epsilon(void){
float fepsilon=1;
double depsilon=1;
long double ldepsilon=1;


while(1+fepsilon!=1){fepsilon/=2;}
fepsilon*=2;
printf("My float epsilon: %g. Library value: %f.\n",fepsilon,FLT_EPSILON);

for (depsilon=1;1+depsilon!=1;depsilon/=2){};
depsilon*=2;
printf("My double epsilon: %g. Library value: %g.\n",depsilon,DBL_EPSILON);

do{ldepsilon/=2;}while(ldepsilon+1!=1);
ldepsilon*=2;
printf("My long double epsilon: %Lg. Library value: %Lg.\n",ldepsilon, LDBL_EPSILON);

}
