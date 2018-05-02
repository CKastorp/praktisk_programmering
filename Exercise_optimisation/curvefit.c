#include"stdio.h"
#include"math.h"
#include"gsl/gsl_multimin.h"
#include"gsl/gsl_vector.h"
#include"gsl/gsl_errno.h"

struct data {int n;double *t, *y, *e;};

double fit(gsl_vector* v,void* params){
double A=gsl_vector_get(v,0), T=gsl_vector_get(v,1), B=gsl_vector_get(v,2);
double S=0;
int i;
struct data *p=(struct data*)params;
int n=p->n;
double* t=p->t;
double* y=p->y;
double* e=p->e;
double err;
for(i=0;i<n;i++){
err=A*exp(-t[i]/T)+B-y[i];
S+=err*err;
}
return S;
}



int curvefit(void){
double t[]={0.47,1.41,2.46,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
double y[]={5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
double e[]={0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
int n=sizeof(t)/sizeof(t[0]);

int d=3, iter=0, status, iteration_status;
double acc=0.01;

struct data p;
p.t=t;
p.y=y;
p.e=e;
p.n=n;

gsl_multimin_fminimizer* system=gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,d);
gsl_multimin_function F;
F.f=fit;
F.n=d;
F.params=(void*)&p;

gsl_vector* start=gsl_vector_alloc(d);
gsl_vector_set(start,0,7);
gsl_vector_set(start,1,1);
gsl_vector_set(start,2,0);

gsl_vector* stepsize=gsl_vector_alloc(d);
gsl_vector_set(stepsize,0,0.2);
gsl_vector_set(stepsize,1,0.2);
gsl_vector_set(stepsize,2,0.2);

gsl_multimin_fminimizer_set(system,&F,start,stepsize);

while(iter<1000){
iter++;
iteration_status=gsl_multimin_fminimizer_iterate(system);
if (iteration_status!=0){fprintf(stderr,"Unable to improve\n");break;}
status=gsl_multimin_test_size(system->size,acc);
if(status ==GSL_SUCCESS){fprintf(stderr,"Converged\n"); break;}
}

fprintf(stderr,"Printing data points\n");
int i;
double A=gsl_vector_get(system->x,0), T=gsl_vector_get(system->x,1), B=gsl_vector_get(system->x,2);
printf("\n\n");
for(i=0;i<n;i++){printf("%g %g \n",t[i],y[i]);}
printf("\n\n");
fprintf(stderr,"Printing fit\n");

double k,f;
for(k=0;k<10;k+=1.0/16){
/*fprintf(stderr,"%g\n",k);*/
f=A*exp(-k/T)+B;
printf("%g %g \n",k,f);}
fprintf(stderr,"print over\n");

gsl_vector_free(stepsize);
gsl_vector_free(start);
gsl_multimin_fminimizer_free(system);

return GSL_SUCCESS;
}
