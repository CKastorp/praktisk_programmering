#include"gsl/gsl_vector.h"

#ifndef HAVE_NEURONS
#define HAVE_NEURONS
typedef struct{int n; double (*f)(double);gsl_vector* data;} ann;
ann* ann_alloc(int n,double f(double));
double ann_feed_forward(ann* network,double x);
void ann_train(ann* network,gsl_vector* x,gsl_vector* y);
void ann_free(ann* network);
int qnewton(double func(gsl_vector*x),gsl_vector* x0,gsl_vector* df,double epsilon);
#endif