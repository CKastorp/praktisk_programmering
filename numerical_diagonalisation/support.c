#include"gsl/gsl_matrix.h"
#include"stdlib.h"

void make_sym_matrix(gsl_matrix* A){
if (A->size1==A->size2){
int n=A->size1,i,j;
for(i=0;i<n;i++){
    gsl_matrix_set(A,i,i,(double) drand48());
    for(j=0;j<i;j++){
    double Aij=(double) drand48();
    gsl_matrix_set(A,i,j,Aij);gsl_matrix_set(A,j,i,Aij);
}}
}}

void print_vector(gsl_vector* v){
for(int i=0;i<v->size;i++) printf("%g ", gsl_vector_get(v,i));
printf("\n\n");
}

void print_matrix(gsl_matrix* Q){
int i,j;
for (i=0;i<Q->size1;i++){
for (j=0;j<Q->size2;j++){
printf("%g ",gsl_matrix_get(Q,i,j));
}
printf("\n");
}
printf("\n");
}

double dotproduct(gsl_vector* A,gsl_vector* B){double result=0;
    for(int i=0;i<A->size;i++){
    double Ai=gsl_vector_get(A,i);
    double Bi=gsl_vector_get(B,i);
result+=Ai*Bi;
}
return result;
}