#include"stdlib.h"
#include"stdio.h"
#include"Nvector.h"

nvector* nvector_alloc(int n){
nvector* v=malloc(sizeof(nvector));
v->size=n;
v->data=malloc(n*sizeof(double));

return v;
}
void nvector_free(nvector* v){
free(v->data);
free(v);
}

void nvector_set(nvector* v, int i, double value){
if(i<(*v).size) (*v).data[i]=value;
else fprintf(stderr, "Index larger than vector length\n");
}

double nvector_get(nvector* v, int i){
if(i<(*v).size) return (*v).data[i];
else fprintf(stderr, "Index larger than vector length\n");
}

double nvector_dot_product(nvector* a, nvector*b){
if (a->size == b->size){
double result=0;
for(int i=0;i<(*a).size;i++) result+=nvector_get(a,i)*nvector_get(b,i);
return result;
}
else {fprintf(stderr, "Error: vectors must be the same size\n");return -1;}
}
