#include"gsl/gsl_errno.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_vector.h"
#include"stdio.h"
#include"math.h"
int QRleastsquare(gsl_vector*x,gsl_vector*y,gsl_vector*dy,gsl_vector*c, double fitfunction(double x,int i));
void loadvectors(gsl_vector*x,gsl_vector*y,gsl_vector*dy);
int QRfactor(gsl_matrix* A,gsl_matrix* R);
int QRsolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

double fitfunction(double x,int i){if (i==0)return log(x);if (i==1)return 1; if(i==2)return x;}

int main(){
    int n=9,k=3;
gsl_vector* x=gsl_vector_alloc(n);
gsl_vector* y=gsl_vector_alloc(n);
gsl_vector* dy=gsl_vector_alloc(n);
gsl_vector* c=gsl_vector_alloc(k);


loadvectors(x,y,dy);/*Hardcoded values*/

int status=QRleastsquare(x,y,dy,c,fitfunction);

printf("Coefficients for log(x),1,x by QR: (%g,%g,%g).\n",gsl_vector_get(c,0),gsl_vector_get(c,1),gsl_vector_get(c,2));

gsl_vector_free(x);
gsl_vector_free(y);
gsl_vector_free(dy);
gsl_vector_free(c);
    return 0;
}