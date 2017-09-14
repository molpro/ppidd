
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

/* C datatypes */
#define PPIDD_C_INT        0
#define PPIDD_C_DOUBLE     1
#define PPIDD_C_LONG       2
#define PPIDD_C_LONG_LONG  3
#define PPIDD_C_FLOAT      4

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
    extern int ProcID(void);
    extern void MPIGA_Error(const char *, int);
    extern void mpi_test_status(const char *, int);
    extern int mpiga_type_f2cmpi(int , MPI_Datatype *);
    extern int mpiga_type_c2cmpi(int , MPI_Datatype *, int *);
    extern int mpiga_type_cmpi2c(MPI_Datatype , int *);
#ifdef __cplusplus
}
#endif

#endif
