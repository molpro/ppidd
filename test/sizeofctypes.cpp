#include <cstdio>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
void sizeofctype(double *sizeint, double *sizelog, double *sizedouble, double *sizefloat)
{
    int verbose=0;
    if (verbose) {
       printf(" The size of some C types on this machine are listed as follows:\n");
       printf("                             int:%5d bytes\n", (int)sizeof(int));
       printf("                            long:%5d bytes\n", (int)sizeof(long));
       printf("                       long long:%5d bytes\n", (int)sizeof(long long));
       printf("                           float:%5d bytes\n", (int)sizeof(float));
       printf("                          double:%5d bytes\n", (int)sizeof(double));
       printf("                         int64_t:%5d bytes\n", (int)sizeof(int64_t));
    }
    *sizeint=(double)sizeof(int64_t);
    *sizelog=(double)sizeof(int);
    *sizedouble=(double)sizeof(double);
    *sizefloat=(double)sizeof(float);
}
#ifdef __cplusplus
}
#endif
