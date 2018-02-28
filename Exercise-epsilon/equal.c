#include<math.h>
#include<stdlib.h>
int equal(double a, double b, double tau, double epsilon){
double absolute = abs(a-b);
double relative = 2*abs(a-b)/(abs(a)+abs(b));
if (absolute < tau){return 1;}
else if (relative < epsilon){
return 1;}
else return 0;
}
