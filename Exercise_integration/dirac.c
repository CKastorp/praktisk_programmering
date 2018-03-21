#include"gsl/gsl_integration.h"
#include"math.h"
#include"stdio.h"
#include"gsl/gsl_errno.h"
double integrand_hamilton(double x, void* params){
double alpha=*(double*)params;
double result =exp(-alpha*x*x)*(-alpha*alpha*x*x/2+alpha/2+x*x/2);
return result;
}
double integrand_norm(double x,void* params){
double alpha=*(double*)params;
return exp(-alpha*x*x);
}
double dirac(double alpha){
double energy;
int limit=1000;

gsl_function f;
f.function=integrand_hamilton;
f.params=(void*)&alpha;
double result_hamilton, result_norm=1, abs=1e-6, rel=1e-6,err;

gsl_integration_workspace* hamiltonian=gsl_integration_workspace_alloc(limit);
int status_1=gsl_integration_qagiu(&f,0,abs,rel,limit,hamiltonian,&result_hamilton,&err);
gsl_integration_workspace_free(hamiltonian);

if (status_1 != GSL_SUCCESS) {return NAN;fprintf(stderr,"Error in hamiltonian integral\n");}
else{
gsl_function g;
g.function=integrand_norm;
g.params=(void*)&alpha;

gsl_integration_workspace* norm=gsl_integration_workspace_alloc(limit);
int status_2=gsl_integration_qagiu(&g,0,abs,rel,limit,norm,&result_norm,&err);
gsl_integration_workspace_free(norm);

if (status_2 != GSL_SUCCESS){return NAN; fprintf(stderr, "Error in norm integral\n");}
else {energy = result_hamilton/result_norm;
return energy;
}
}
}
