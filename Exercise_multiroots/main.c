#include"stdio.h"
#include"gsl/gsl_vector.h"
#include"stdlib.h"
void findmin(gsl_vector*);
int hydrogen_master(double);
int main(){
gsl_vector* result=gsl_vector_alloc(2);

findmin(result);
printf("Rosenbrock minimum: (%g,%g).\n",gsl_vector_get(result,0),gsl_vector_get(result,1));
printf("Analytic result: (1,1)\n");

double rmax=8;
int flag=hydrogen_master(rmax);

printf("\n\n");
rmax=10;
flag=hydrogen_master(rmax);
gsl_vector_free(result);


if (flag==GSL_SUCCESS) return 0;
else{fprintf(stderr,"Error: %s\n",gsl_strerror(flag)); return flag;}
}
