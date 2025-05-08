/*! PPIDD standalone c test suites */
#include "ppidd.h"
#include <stdio.h>
#include <inttypes.h>

int ppidd_ctest();

int main(int argc, char **argv)
{
    /* Initialize PPIDD */
    PPIDD_Initialize(&argc, &argv, PPIDD_IMPL_DEFAULT, 32);
    PPIDD_Initialize_data();

    ppidd_ctest();

    /* Terminate and Tidy up PPIDD */
    PPIDD_Finalize();
    return 0;
}

int ppidd_ctest()
{
    int nproc = PPIDD_Size();
    int me = PPIDD_Rank();

    if(me==0) {
       printf(" PPIDD initialized\n");
       printf(" Nprocs= %d    My proc= %d\n",nproc,me);
       printf(" Performing PPIDD C tests:\n");
       fflush(stdout);
    }

    PPIDD_Barrier();

    if(me==0) {
       printf(" All PPIDD C tests successful.\n");
       fflush(stdout);
    }

    return 0;
}
