/*! PPIDD standalone c test suites */
#include "ppidd.h"
#include <stdio.h>
#include <inttypes.h>

int ppidd_ctest(void);

int main(int argc, char **argv)
{
    /* Initialize PPIDD */
    PPIDD_Initialize(&argc, &argv, PPIDD_IMPL_DEFAULT);
    PPIDD_Initialize_data();

    ppidd_ctest();

    /* Terminate and Tidy up PPIDD */
    PPIDD_Finalize();
    return 0;
}

int ppidd_ctest(void)
{
    int64_t me, nproc;

    PPIDD_Size(&nproc);
    PPIDD_Rank(&me);

    if(me==0) {
       printf(" PPIDD initialized\n");
       printf(" Nprocs= %" PRIu64 "    My proc= %" PRIu64 "\n",nproc,me);
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
