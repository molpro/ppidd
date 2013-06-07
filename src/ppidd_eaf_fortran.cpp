
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
 * FORTRAN interface.  The subroutines in this file named PPIDD_XXXXX are *
 * converted to the proper FORTRAN external by the FORT_Extern macro and  *
 * the definitions in the ppidd_eaf_fortran.h header file.                *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       15/07/2008                                                 *
\* ====================================================================== */


#include "ppidd_eaf_fortran.h"   /* include ppidd_machines.h */

#ifdef MPI2
 #include <mpi.h>
 #define   MPI_EAF_RW -1
 #define   MPI_EAF_W  -2
 #define   MPI_EAF_R  -3
 extern MPI_Comm MPIGA_WORK_COMM;
#endif

#ifdef GA_MPI
 #include <ga.h>
 #include <eaf.h>
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
