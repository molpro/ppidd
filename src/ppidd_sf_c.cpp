
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
 * C interface.  The subroutines in this file named PPIDD_XXXXX can be    *
 * only called by C program directly.  Any calling by Fortran progam      *
 * should refer to the routines in the ppidd_sf_fortran.h header file.    *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       15/07/2008                                                 *
\* ====================================================================== */

#ifdef MPI2
 #include "ppidd_mpi.h"
 extern MPI_Comm MPIGA_WORK_COMM;
#endif

#include "ppidd_sf_c.h"   /* include ppidd_machines.h */

#ifdef FORTCL_NEXT
#undef FORTCL_NEXT
#endif

#ifdef FORTCL_END
#undef FORTCL_END
#endif

#ifdef FORTINTC_DIVIDE
#undef FORTINTC_DIVIDE
#endif

/* The following code should be the same as those in ppidd_sf_fortran.cpp (except ppidd_sf_rank). One should make it consistent once code in ppidd_sf_fortran.cpp is changed. */



#ifdef GA_MPI
 #include <ga.h>
 #include <sf.h>
#endif
 extern int ppidd_sf_rank(void);
#if defined(GA_MPI) || defined(MPI2)
 static int MPI_Debug=0;
#endif

/* Common code shared by Fortran and C interfaces for PPIDD Shared Files Library.*/
#include "ppidd_sf_share.h"
