#include"math.h"
#include"gsl/gsl_integration.h"
#include"gsl/gsl_errno.h"
double integrand(double x, void* params){
return log(x)/sqrt(x);

}
double integration(){
double result,acc=1e-6,rel=1e-6,err;

gsl_function f;
f.function=integrand;

int limit=100;
gsl_integration_workspace* workspace=gsl_integration_workspace_alloc(limit);
int flag=gsl_integration_qags(&f,0,1,acc,rel,limit,workspace,&result,&err);
gsl_integration_workspace_free(workspace);
if (flag !=GSL_SUCCESS) return NAN;
else return result;
}
