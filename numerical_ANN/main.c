#include "ANN.h"
#include"gsl/gsl_errno.h"
#include"stdio.h"
#include"math.h"
#include"gsl/gsl_vector.h"

double gauss(double x) {return exp(-x*x);}
double testfunction(double x) {return sin(x)*exp(-x*x);}
double secondtest(double x){return cos(5*x)*exp(-fabs(x));}
int main(){
    
int n=4;
ann* network=ann_alloc(n,gauss);
double a=-5,b=5;
int nx=100;
gsl_vector* vx=gsl_vector_alloc(nx);
gsl_vector* vf=gsl_vector_alloc(nx);

for(int i=0;i<nx;i++){
    double x=a+(b-a)*i/(nx-1);
    gsl_vector_set(vx,i,x);
    gsl_vector_set(vf,i,testfunction(x));
}
for(int i=0;i<n;i++){
gsl_vector_set(network->data,i*3,a+(b-a)*i/(n-1));
gsl_vector_set(network->data,3*i+1,1);
gsl_vector_set(network->data,3*i+2,pow(-1,i));
}

ann_train(network,vx,vf);

for(int i=0;i<nx;i++)printf("%g %g\n",gsl_vector_get(vx,i),gsl_vector_get(vf,i));
printf("\n\n");
for(double k=a;k<b;k+=1.0/16){printf("%g %g\n",k,ann_feed_forward(network,k));
}

ann_free(network);
gsl_vector_free(vf);
gsl_vector_free(vx);
    return 0;
}