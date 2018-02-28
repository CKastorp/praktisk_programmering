#include "stdio.h"
#include "math.h"
#include "limits.h"
#include "float.h"
void maxima(){
int i=0;

while (i+1>i){i++;}
printf("My maximum integer is %i.\n",i);

i=0;

while (i-1<i){i=i-1;}
printf("My smallest integer is %i.\n",i);

for (i=0;i+1>i;i++){}
printf("Largest integer, for loop: %i\n",i);

for (i=0;i-1<i;i++){}
printf("Smallest integer, for loop: %i\n",i);

i=0;
do{i++;}while(i+1>i);
printf("Largest integer do-while: %i\n",i);

i=0;
do{i++;}while(i-1<i);
printf("Smallest integer, do-while: %i\n",i);

}
