#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

/* ====================================================================== *\
 *                 PPIDD Exclusive Access File Library                    *
 *                 ===================================                    *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ------------------------------------------------------------------------------------------ *
 * An exclusive access file is a file which is generated and/or read by a single process of a *
 * distributed parallel application. Files are not shared between different processes. The    *
 * library is an abstract high-performance file system which provides a common interface for  *
 * a variety of architecture specific parallel storage systems.  The library also makes       *
 * available features like asynchronous input and output to Fortran.  EAF's syntax is similar *
 * to the standard Unix C file operations, differences indicate new semantics or extended     *
 * features available through EAF.                                                            *
 *    The last argument of all subroutines returns an integer error code with the value zero  *
 * implying success, non-zero implying some error condition.  Offsets are doubles and an      *
 * offset with a fractional component generates an error.                                     *
 *--------------------------------------------------------------------------------------------*
 * FORTRAN interface.
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       15/07/2008                                                 *
\* ====================================================================== */




#ifdef MPI2
 #include <mpi.h>
 #define   MPI_EAF_RW -1
 #define   MPI_EAF_W  -2
 #define   MPI_EAF_R  -3
 extern MPI_Comm MPIGA_WORK_COMM;
#endif

/* Map the function names to something the Fortran compiler can call */
#define PPIDD_Eaf_open     FC_FUNC_(ppidd_eaf_open,PPIDD_EAF_OPEN)
#define PPIDD_Eaf_write    FC_FUNC_(ppidd_eaf_write,PPIDD_EAF_WRITE)
#define PPIDD_Eaf_awrite   FC_FUNC_(ppidd_eaf_awrite,PPIDD_EAF_AWRITE)
#define PPIDD_Eaf_read     FC_FUNC_(ppidd_eaf_read,PPIDD_EAF_READ)
#define PPIDD_Eaf_aread    FC_FUNC_(ppidd_eaf_aread,PPIDD_EAF_AREAD)
#define PPIDD_Eaf_wait     FC_FUNC_(ppidd_eaf_wait,PPIDD_EAF_WAIT)
#define PPIDD_Eaf_waitall  FC_FUNC_(ppidd_eaf_waitall,PPIDD_EAF_WAITALL)
#define PPIDD_Eaf_probe    FC_FUNC_(ppidd_eaf_probe,PPIDD_EAF_PROBE)
#define PPIDD_Eaf_close    FC_FUNC_(ppidd_eaf_close,PPIDD_EAF_CLOSE)
#define PPIDD_Eaf_delete   FC_FUNC_(ppidd_eaf_delete,PPIDD_EAF_DELETE)
#define PPIDD_Eaf_length   FC_FUNC_(ppidd_eaf_length,PPIDD_EAF_LENGTH)
#define PPIDD_Eaf_truncate FC_FUNC_(ppidd_eaf_truncate,PPIDD_EAF_TRUNCATE)
#define PPIDD_Eaf_errmsg   FC_FUNC_(ppidd_eaf_errmsg,PPIDD_EAF_ERRMSG)

#ifdef GA_MPI
 #include <ga.h>
#ifdef __cplusplus
extern "C" {
#endif
 #include <eaf.h>
#ifdef __cplusplus
}
#endif
#endif

#if defined(GA_MPI) || defined(MPI2)
 static int MPI_Debug=0;
#endif


/* ************************************************************************ *\
   Get calling process id in the work communicator. This function is only
   used in operations related to eaf.
\* ************************************************************************ */
   int ppidd_eaf_rank(void) {
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

/* Common code shared by Fortran and C interfaces for PPIDD EAF Library.*/
#ifdef __cplusplus
extern "C" {
#endif
#include "ppidd_eaf_share.h"
#ifdef __cplusplus
}
#endif
