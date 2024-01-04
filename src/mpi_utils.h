
/* ==================================================================== *\
 * Standalone Utility Tools for MPI
\* ==================================================================== */

/* ------------- *\
***Include Files***
\* ------------- */

#ifndef __MPI_UTILS_H__
#define __MPI_UTILS_H__

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* some definitions for mpi_test_status */
#if !defined PPIDD_MIN
  #define PPIDD_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#endif
#define MSG_ERR_STR_LEN 80
#define TEST_ERR_STR_LEN (MSG_ERR_STR_LEN+MPI_MAX_ERROR_STRING)
extern  char  mpi_test_err_string[TEST_ERR_STR_LEN];



/* =========================== *\
    mpi_utils Function Prototypes
\* =========================== */

#ifdef __cplusplus
extern "C" {
#endif
    extern int NNodes_Total(MPI_Comm, int *);
    extern int ProcID();
    extern void MPIGA_Error(const char *, int);
    extern void mpi_test_status(const char *, int);
    extern MPI_Datatype dtype_mpi(int);
    extern size_t dtype_size(int);
#ifdef __cplusplus
}
#endif

#endif
