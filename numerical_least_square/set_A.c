#include"gsl/gsl_matrix.h"
#include"math.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_errno.h"

int hardcode_fitfunction(gsl_vector* x,gsl_matrix* A){
if(x->size !=A->size1) return GSL_EINVAL;
else{
int i,n=x->size;
double xi;
for(i=0;i<n;i++){
xi=gsl_vector_get(x,i);
gsl_matrix_set(A,i,0,log(xi));
gsl_matrix_set(A,i,1,1);
gsl_matrix_set(A,i,2,xi);
}

return GSL_SUCCESS;
}}

void loadvectors(gsl_vector* x,gsl_vector* y,gsl_vector* dy){
    gsl_vector_set(x,0,0.1); gsl_vector_set(x,1,1.33);gsl_vector_set(x,2,2.55);
    gsl_vector_set(x,3,3.78);gsl_vector_set(x,4,5);gsl_vector_set(x,5,6.22);
    gsl_vector_set(x,6,7.45);gsl_vector_set(x,7,8.68);gsl_vector_set(x,8,9.9);

    gsl_vector_set(y,0,-15.3);gsl_vector_set(y,1,0.32);gsl_vector_set(y,2,2.45);
    gsl_vector_set(y,3,2.75);gsl_vector_set(y,4,2.27);gsl_vector_set(y,5,1.35);
    gsl_vector_set(y,6,0.157);gsl_vector_set(y,7,-1.23);gsl_vector_set(y,8,-2.75);

    gsl_vector_set(dy,0,1.04);gsl_vector_set(dy,1,0.594);gsl_vector_set(dy,2,0.983);
    gsl_vector_set(dy,3,0.998);gsl_vector_set(dy,4,1.11);gsl_vector_set(dy,5,0.398);
    gsl_vector_set(dy,6,0.535);gsl_vector_set(dy,7,0.968);gsl_vector_set(dy,8,0.478);
}