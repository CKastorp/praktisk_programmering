#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_sf_airy.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
double x;
void airy(void);
int main(){
airy();
double matrix_a[]={6.13, -2.9, 5,80,
		8.08, -6.31, -3.89,
		-4.36, 1.00, 0.19};
double safe_copy[]={6.13, -2.9, 5,80,
                8.08, -6.31, -3.89,
                -4.36, 1.00, 0.19};
double RightHandSide[]={6.23, 5.37, 2.29};
int sign;
gsl_matrix_view M=gsl_matrix_view_array(matrix_a,3,3);
gsl_vector_view b=gsl_vector_view_array(RightHandSide,3);
gsl_vector* x=gsl_vector_alloc(3);
gsl_permutation* P=gsl_permutation_alloc(3);

gsl_linalg_LU_decomp(&M.matrix,P,(&sign));
gsl_linalg_LU_solve(&M.matrix,P,&b.vector,x);
gsl_permutation_free(P);

printf("Solution 'x' from LU_solve: (%g,%g,%g)\n",gsl_vector_get(x,0)\
	,gsl_vector_get(x,1),gsl_vector_get(x,2));

double X= safe_copy[0]*gsl_vector_get(x,0)+safe_copy[1]*gsl_vector_get(x,1)\
	+safe_copy[2]*gsl_vector_get(x,2);
double Y= safe_copy[3]*gsl_vector_get(x,0)+safe_copy[4]*gsl_vector_get(x,1)\
	+safe_copy[5]*gsl_vector_get(x,2);
double Z= safe_copy[6]*gsl_vector_get(x,0)+safe_copy[7]*gsl_vector_get(x,1)\
	+safe_copy[8]*gsl_vector_get(x,2);
printf("Solution times coefficient matrix: (%g,%g,%g)\n",X,Y,Z);

printf("Original right hand side: (%g,%g,%g)\n",RightHandSide[0],\
	RightHandSide[1],RightHandSide[2]);
gsl_vector_free(x);
return 0;}
