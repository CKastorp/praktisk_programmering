#include"komplex.h"
#include"stdio.h"
#include"math.h"
void komplex_print(char *s, komplex a){
printf("%s %g%+gi\n", s,a.re, a.im);}

komplex komplex_new (double x, double y){
komplex z={x,y};
return (z);
}

void komplex_set(komplex* z, double x, double y){
z->re=x;
z->im=y;
}

komplex komplex_add(komplex a,komplex b){
komplex result={a.re+b.re,a.im+b.im};
return(result);
}

int komplex_equal(komplex a,komplex b){
int i=0;
if (a.re==b.re){
	if (a.im==b.im){i=1;}
}
return(i);
}


komplex komplex_mul(komplex a, komplex b){
double real = a.re*b.re-a.im*b.im;
double imag = a.im*b.re+a.re*b.im;
komplex result ={real,imag};
return(result);
}

komplex komplex_conjugate(komplex a){
komplex result={a.re,-a.im};
return(result);}

komplex komplex_div(komplex a, komplex b){
komplex nominator = komplex_mul(a,komplex_conjugate(b));
double denominator = pow(b.re,2)+pow(b.im,2);
komplex output={nominator.re/denominator,nominator.im/denominator};
return(output);
}
