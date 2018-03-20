#ifndef HAVE_NVECTOR_H

typedef struct {int size; double* data;} nvector;

nvector* nvector_alloc       (int n);                           /* allocates memory for size-n vector */
void     nvector_free        (nvector* v);                      /* frees memory */
void     nvector_set         (nvector* v, int i, double value); /* v_{i} = value */;
double   nvector_get         (nvector* v, int i);               /* returns v_{i} */
double   nvector_dot_product (nvector* u, nvector* v);          /* dot-product */

/* optional */
void nvector_print    (char* s, nvector* v);    /* prints s and then the vector */
int  double_equal     (double a, double b);     /* 1 if equal, 0 otherwise */
int  nvector_equal    (nvector* a, nvector* b);     /* 1 if equal, 0 if not */
void nvector_set_zero (nvector* v);             /* all elements = 0 */
void nvector_add      (nvector* a, nvector* b); /* a_i <- a_i + b_i */
void nvector_sub      (nvector* a, nvector* b); /* a_i <- a_i - b_i */
void nvector_scale    (nvector* a, double x);   /* a_i <- x*a_i     */


#define HAVE_NVECTOR_H
#endif
