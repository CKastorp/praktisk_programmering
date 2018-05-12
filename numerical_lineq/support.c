#include"gsl/gsl_matrix.h"
#include"stdio.h"
#include"gsl/gsl_vector.h"
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

double matrix_column_dotproduct(gsl_matrix* A,int i,gsl_matrix* B,int j){
double result=0;
for(int n=0;n<A->size1;n++)result+=gsl_matrix_get(A,n,i)*gsl_matrix_get(B,n,j);

    return result;
}

void matrixproduct(gsl_matrix* A,gsl_matrix* B, gsl_matrix* result){/*matrix product A*B*/
double temp;
for (int i=0;i<B->size2;i++){
for (int j=0;j<A->size1;j++){
    temp=0;
    for (int k=0;k<B->size1;k++)temp+=gsl_matrix_get(B,k,i)*gsl_matrix_get(A,j,k);
    gsl_matrix_set(result,j,i,temp);
}
}
}
void matrixtranspose(gsl_matrix* A, gsl_matrix* AT){
for(int i=0;i<A->size1;i++){
    for(int j=0;j<A->size2;j++){
        double x=gsl_matrix_get(A,i,j);
        gsl_matrix_set(AT,j,i,x);}
}

}

void random_matrix(gsl_matrix* A){
for(int i=0;i<A->size1;i++){
for(int j=0;j<A->size2;j++)
gsl_matrix_set(A,i,j,(double) rand()/RAND_MAX*20);
}

}