#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"gsl/gsl_errno.h"
#include"math.h"
double dotproduct(gsl_vector* A,gsl_vector* B);
double matrix_column_dotproduct(gsl_matrix* A, int i, gsl_matrix* B, int j);

int QRfactor(gsl_matrix* A,gsl_matrix* R){
int m=A->size2;
int n=A->size1;

if (R->size1!=R->size2)return GSL_ENOTSQR;
if (R->size1!=m)return GSL_EINVAL;

gsl_matrix_set_zero(R);

double overlap;
double norm;
double dummy;

for (int i=0;i<m;i++){
norm=sqrt(matrix_column_dotproduct(A,i,A,i));
gsl_matrix_set(R,i,i,norm);
for(int j=0;j<n;j++){dummy=gsl_matrix_get(A,j,i)/norm;
gsl_matrix_set(A,j,i,dummy);}

for(int k=i+1;k<m;k++){
dummy=matrix_column_dotproduct(A,i,A,k);
gsl_matrix_set(R,i,k,dummy);

for(int l=0;l<n;l++){
    overlap=gsl_matrix_get(A,l,k)-dummy*gsl_matrix_get(A,l,i);
    gsl_matrix_set(A,l,k,overlap);}
}
}

    return GSL_SUCCESS;
}