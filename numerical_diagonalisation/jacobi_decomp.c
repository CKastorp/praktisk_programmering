#include"stdio.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_errno.h"
#include"gsl/gsl_blas.h"
#include"math.h"
double dotproduct(gsl_vector* A,gsl_vector* B);

int eig_sweep(gsl_matrix* A,gsl_vector* values, gsl_matrix* V){
if (A->size1!=A->size2) return GSL_ENOTSQR;
if (V->size1!=V->size2) return GSL_ENOTSQR;
if (A->size1!=V->size1) return GSL_EBADLEN;
if (A->size1!=values->size) return GSL_EBADLEN;

int n=A->size1;
gsl_matrix_set_identity(V);
for(int l=0;l<n;l++)gsl_vector_set(values,l,gsl_matrix_get(A,l,l));

int iter=0,change;
/*double err=0.0001,dev=1;*/
while(iter<200){
    iter++;
change=0;

for(int p=0;p<n;p++){
    
for (int q=p+1;q<n;q++){
double App=gsl_vector_get(values,p);
    double Aqq=gsl_vector_get(values,q);
        double Apq=gsl_matrix_get(A,p,q);
double phi=0.5*atan2(2*Apq,Aqq-App);
double c=cos(phi),s=sin(phi);

for(int i=0;i<n;i++) {double Vip=gsl_matrix_get(V,i,p),Viq=gsl_matrix_get(V,i,q);
gsl_matrix_set(V,i,p,c*Vip-s*Viq);gsl_matrix_set(V,i,q,s*Vip+c*Viq);}
double App_prime=App*c*c+s*s*Aqq-2*s*c*Apq;
double Aqq_prime=s*s*App+Aqq*c*c+2*s*c*Apq;
if(App_prime !=App) change=1;
if (Aqq_prime !=Aqq) change=1;
/*double Apq_prime=s*c*(App-Apq)+(c*c-s*s)*Apq;*/
gsl_vector_set(values,p,App_prime);gsl_vector_set(values,q,Aqq_prime);
gsl_matrix_set(A,p,q,0.0);
for (int i=0;i<p;i++){
double Aip=gsl_matrix_get(A,i,p);double Aiq=gsl_matrix_get(A,i,q);
gsl_matrix_set(A,i,p,c*Aip-s*Aiq);gsl_matrix_set(A,i,q,s*Aip+c*Aiq);}
for(int i=p+1;i<q;i++){
double Api=gsl_matrix_get(A,p,i);double Aiq=gsl_matrix_get(A,i,q);
gsl_matrix_set(A,p,i,c*Api-s*Aiq);gsl_matrix_set(A,i,q,s*Api+c*Aiq);}

for(int i=q+1;i<n;i++){
double Api=gsl_matrix_get(A,p,i);double Aqi=gsl_matrix_get(A,q,i);
gsl_matrix_set(A,p,i,c*Api-s*Aqi);gsl_matrix_set(A,q,i,s*Api+c*Aqi);}
}}


if(change==0)break;
}

printf("Number of sweeps: %i.\n\n",iter);
/*for (int i=0;i<n;i++){
    for(int j=i+1;j<n;j++) gsl_matrix_set(A,i,j,gsl_matrix_get(A,j,i));
}*/


    return GSL_SUCCESS;
}