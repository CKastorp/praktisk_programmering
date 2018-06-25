#include "ANN.h"
#include"gsl/gsl_errno.h"
#include"stdio.h"
#include"math.h"
#include"gsl/gsl_vector.h"

ann* ann_alloc(int n,double f(double)){
    ann* pointer=malloc(sizeof(ann));
    pointer->n=n;
    pointer->f=f;
    pointer->data=gsl_vector_alloc(3*n);
    return pointer;
}

double ann_feed_forward(ann* network,double x){
    double s=0;
    double n=network->n;
    for(int i=0;i<n;i++){
        double a=gsl_vector_get(network->data,3*i);
        double b=gsl_vector_get(network->data,3*i+1);
        double w=gsl_vector_get(network->data,i*3+2);
        s+=network->f((x-a)/b)*w;}
        return s;
}



void ann_train(ann* network,gsl_vector* x,gsl_vector* y){

double neuron_error(gsl_vector* p){
gsl_vector_memcpy(network->data,p);
double sum=0;
int L=x->size;
for(int i=0;i<L;i++){
double x0=gsl_vector_get(x,i);
double y0=gsl_vector_get(y,i);
double fx=ann_feed_forward(network,x0);
sum+=pow(y0-fx,2);
}
return sum/L;
}
int n=network->data->size;
gsl_vector* p=gsl_vector_alloc(n);
gsl_vector* dp=gsl_vector_alloc(n);
gsl_vector_memcpy(p,network->data);
qnewton(neuron_error,p,dp,1e-5);
gsl_vector_memcpy(network->data,p);

gsl_vector_free(p);
gsl_vector_free(dp);
}
void ann_free(ann* network){

    gsl_vector_free(network->data);
    free(network);
}