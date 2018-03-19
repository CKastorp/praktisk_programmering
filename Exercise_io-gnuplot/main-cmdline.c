#include"math.h"
#include"stdlib.h"
#include"stdio.h"
int main(int argc, char** argv){

if (argc >1){
double x;
int i;
for (i=1;i<argc;i++){
x=atof(argv[i]);
printf("%g \t %g \n", x, sin(x));
}
}
else printf("%i \t %i \n",0,0);
return 0;
}
