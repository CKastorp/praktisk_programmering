#include"Nvector.h"
#include"stdio.h"
#include"stdlib.h"
int main(){
int n=3;
printf("\nTesting nvec_alloc\n");
nvector* a=nvector_alloc(n);
nvector* b=nvector_alloc(n);
printf("Test succesful\n");

printf("\nTesting nvector set/get\n");
double b0=4;
double b1=5;
double b2=1;
for(int i=0;i<n;i++)nvector_set(a,i,i);
nvector_set(b,0,b0);nvector_set(b,1,b1);nvector_set(b,2,b2);
double a1_test=nvector_get(a,1);
double b2_test=nvector_get(b,2);
if (1==a1_test) printf("Value is %lg. Should be 1. Test passed\n", a1_test);
if (b2==b2_test) printf("Value is %lg. Should be %g. Test passed\n", b2_test,b2);
else printf("Value is %g. Should be %g. Test failed. Please debug.\n",b2_test,b2);

printf("\nTesting dot product\n");
double real_value=0*b0+1*b1+2*b2;
double test_dot=nvector_dot_product(a,b);
if (test_dot==real_value) printf("Dot product is %g. Should be %g. Test passed.\n",test_dot,real_value);
else printf("Dot product is %g. Should be %g. Test failed.\n",test_dot,real_value);
nvector_free(b);
nvector_free(a);
return 0;
}
