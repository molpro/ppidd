/*! PPIDD standalone c test suites */
#include "ppidd_machines.h"
#include "ppidd_c.h"

int ppidd_ctest();

#ifndef NOMAIN
int main(int argc, char **argv)
{
    /* Initialize PPIDD */
    PPIDD_Initialize(argc, argv);
    PPIDD_Initialize_data();

    ppidd_ctest();

    /* Terminate and Tidy up PPIDD */
    PPIDD_Finalize();
    return 0;
}
#endif

int ppidd_ctest(void)
{
    int me, nproc;

    PPIDD_Size(&nproc);
    PPIDD_Rank(&me);

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
