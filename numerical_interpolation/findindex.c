#include"stdio.h"
#include"math.h"

int binsearch(double x[], double z,int n){
int lower=0;
int upper=n-1;
int index;
if (n<=0){fprintf(stderr, "Error: n<=0\n");return 0;}
if (z<x[0]){fprintf(stderr,"Error: z<x[0]\n"); return 0;}

while (upper-lower>1){index=lower+(upper-lower)/2;
if(z==x[index]) {lower=index;break;}
if(z<x[index])upper=index;
else lower=index;

}
return lower;
}
