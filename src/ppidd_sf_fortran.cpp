#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

/* ====================================================================== *\
 *                    PPIDD Shared Files Library                          *
 *                    ==========================                          *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ------------------------------------------------------------------------------------------ *
 * The Shared Files (SF) library implements logically-shared temporary files for parallel     *
 * SPMD (single-program-multiple-data) applications. The main features are listed as follows: *
 * -- Shared files are non-persistent (temporary)                                             *
 * -- Shared files resemble one-dimensional arrays in main memory                             *
 * -- Each process can independently read/write to any location in the file                   *
 * -- All routines return error code: "0" means success.                                      *
 * -- sf_create and sf_destroy are collective                                                 *
 * -- file, request sizes, and offset (all in bytes) are DOUBLE PRECISION arguments,          *
 *    all the other arguments are INTEGERS                                                    *
 * -- read/writes are asynchronous                                                            *
 *--------------------------------------------------------------------------------------------*
 * FORTRAN interface.  The subroutines in this file named PPIDD_XXXXX are *
 * converted to the proper FORTRAN external by the FC_FUNC macro.         *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       15/07/2008                                                 *
\* ====================================================================== */




#ifdef MPI2
 #include <mpi.h>
 extern MPI_Comm MPIGA_WORK_COMM;
#endif

/* Map the function names to something the Fortran compiler can call */
#define PPIDD_Sf_create  FC_FUNC_(ppidd_sf_create,PPIDD_SF_CREATE)
#define PPIDD_Sf_write   FC_FUNC_(ppidd_sf_write,PPIDD_SF_WRITE)
#define PPIDD_Sf_read    FC_FUNC_(ppidd_sf_read,PPIDD_SF_READ)
#define PPIDD_Sf_wait    FC_FUNC_(ppidd_sf_wait,PPIDD_SF_WAIT)
#define PPIDD_Sf_waitall FC_FUNC_(ppidd_sf_waitall,PPIDD_SF_WAITALL)
#define PPIDD_Sf_destroy FC_FUNC_(ppidd_sf_destroy,PPIDD_SF_DESTROY)
#define PPIDD_Sf_errmsg  FC_FUNC_(ppidd_sf_errmsg,PPIDD_SF_ERRMSG)

#ifdef GA_MPI
 #include <ga.h>
#ifdef __cplusplus
extern "C" {
#endif
 #include <sf.h>
#ifdef __cplusplus
}
#endif
#endif

#if defined(GA_MPI) || defined(MPI2)
 static int MPI_Debug=0;
#endif


/* ************************************************************************ *\
   Get calling process id in the work communicator. This function is only
   used in operations related to sf.
\* ************************************************************************ */
   int ppidd_sf_rank(void) {
      int myid=0;
#ifdef MPI2
      MPI_Comm mpicomm=MPIGA_WORK_COMM;

      MPI_Comm_rank(mpicomm,&myid);
#endif
#ifdef GA_MPI
      myid=GA_Nodeid();
#endif
      return(myid);
   }

/* Common code shared by Fortran and C interfaces for PPIDD Shared Files Library.*/
#ifdef __cplusplus
extern "C" {
#endif
#include "ppidd_sf_share.h"
#ifdef __cplusplus
}
#endif
