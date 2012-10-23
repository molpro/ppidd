
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
 * converted to the proper FORTRAN external by the FORT_Extern macro and  *
 * the definitions in the ppidd_sf_fortran.h header file.                 *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       15/07/2008                                                 *
\* ====================================================================== */


#include "ppidd_sf_fortran.h"   /* include ppidd_machines.h */

#ifdef MPI2
 #include <mpi.h>
 extern MPI_Comm MPIGA_WORK_COMM;
#endif

#ifdef GA_TOOLS
 #include <ga.h>
 #include <sf.h>

 #ifndef GA_VERSION_GE_5
 extern Integer sf_create(char *fname, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle);
 extern Integer FATR sf_write_(Integer *s_a, SFsize_t *offset, SFsize_t *bytes, char *buffer, Integer *req_id);
 extern Integer FATR sf_read_(Integer *s_a, SFsize_t *offset, SFsize_t *bytes, char *buffer, Integer *req_id);
 extern Integer FATR sf_wait_(Integer *req_id);
 extern Integer FATR sf_waitall_(Integer *list, Integer *num);
 extern Integer FATR sf_destroy_(Integer *s_a);
 extern void sf_errmsg(int code, char *msg);
 #endif

#endif

#if defined(GA_TOOLS) || defined(MPI2)
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
#ifdef GA_TOOLS
      myid=GA_Nodeid();
#endif
      return(myid);
   }

/* Common code shared by Fortran and C interfaces for PPIDD Shared Files Library.*/
#include "ppidd_sf_share.h"
