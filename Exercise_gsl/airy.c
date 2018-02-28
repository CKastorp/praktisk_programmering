#include<stdio.h>
#include<gsl/gsl_sf_airy.h>
#include<tgmath.h>
void airy(){
double x;
double A;
double B;
FILE* destination=fopen("airy.data","w");

gsl_mode_t mode=GSL_PREC_DOUBLE;
for(x=-10;x<10;x+=0.1){
A=gsl_sf_airy_Ai(x,mode);
B=gsl_sf_airy_Bi(x,mode);
fprintf(destination, "%g\t  %g\t  %g \n",x,A,B);}
fclose(destination);

}
