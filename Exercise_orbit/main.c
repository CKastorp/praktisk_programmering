#include"stdio.h"
#include"gsl/gsl_errno.h"

int logistic(void);
int orbit (void);
int main(){
int status;
status = logistic();
if(status != GSL_SUCCESS){printf("Error, return value %i.\n",status);return 1;}
else{printf("\n\n");
status= orbit();
if(status != GSL_SUCCESS){printf("Error, return value %i.\n",status);return 1;}
else{
return 0;
}}
}
