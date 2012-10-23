
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

/* The following code should be the same as those in ppidd_sf_fortran.c (except ppidd_sf_rank). One should make it consistent once code in ppidd_sf_fortran.c is changed. */

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
 extern int ppidd_sf_rank(void);
#if defined(GA_TOOLS) || defined(MPI2)
 static int MPI_Debug=0;
#endif

/* Common code shared by Fortran and C interfaces for PPIDD Shared Files Library.*/
#include "ppidd_sf_share.h"
