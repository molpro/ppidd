#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

/* ====================================================================== *\
 *                 MPI Version of Mutual Exclusion(Mutex)                 *
 *                 ======================================                 *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ---------------------------------------------------------------------- *
 * C sorce code of MPI Version of Mutual exclusion. The subroutines can be*
 * called directly by C code, while the corresponding Fortran wrappers    *
 * (which can be called by Fortran code) are in other files.              *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       30/01/2009                                                 *
\* ====================================================================== */

#ifdef HAVE_MPI_H

#include "mpi_utils.h"
#include "mpiga_base.h"
#include "mpi_nxtval.h"

mpimutex_t_index *onesided_helpmutex_data_struc=NULL, *onesided_helpmutex_index=NULL;
int onesided_helpmutex_num=0;

int use_onesided_helpmutex=0;

static int MPI_Debug = 0;

/* According the usage range, there are two kind of mutexes: one is used in each global data structure; another is used in other general code.
MAX_ONESIDED_HELPMUTEX_ARRAYS=MAX_ONESIDED_HELPMUTEX_GA (used in GA) + MAX_ONESIDED_HELPMUTEX_NONGA (used in other code) */
/* There are two storage methods for mutex window: (1) in computation process; (2) in helper process.
   Here for the routines in this file the mutexes are stored in the helper process. */

/* initialization for one-sided helpmutexes -- called just before install_twosided_nxtval */
void initialize_onesided_helpmutexes()
{
    int i,mpierr;
    int homerank;
    mpimutex_t mutex;
    int  server = LastServerID();         /* id of last server process */

    if (use_helper_server) {
/* zero in pointers in onesided_helpmutex array, initialise onesided_helpmutex data structure */
    if (MPI_Debug) printf("%5d: In initialize_onesided_helpmutexes: sizeof(mpimutex_t_index)=%d, sizeof(mpimutex)=%d, MAX_ONESIDED_HELPMUTEX_ARRAYS=%d\n",
        ProcID(),(int)sizeof(mpimutex_t_index),(int)sizeof(struct mpimutex),MAX_ONESIDED_HELPMUTEX_ARRAYS);
    onesided_helpmutex_data_struc=(mpimutex_t_index *)malloc(sizeof(mpimutex_t_index)*MAX_ONESIDED_HELPMUTEX_ARRAYS);
    if(!onesided_helpmutex_data_struc){
      fprintf(stderr,"ERROR in initialize_onesided_helpmutexes: failed to malloc onesided_helpmutex data structure.\n");
      exit (1);
    }
    onesided_helpmutex_index = onesided_helpmutex_data_struc;
    for(i=0;i<MAX_ONESIDED_HELPMUTEX_ARRAYS; i++) {
       homerank=server;
       mpierr = MPIMUTEX_Create(homerank, MPI_COMM_WORLD, &mutex);
       if (MPI_Debug) printf("%5d: In initialize_onesided_helpmutexes: after calling MPIMUTEX_Create. i=%d\n",ProcID(),i);
       if (mpierr != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
       onesided_helpmutex_index[i].actv = 0;
       onesided_helpmutex_index[i].lock = 0;
       onesided_helpmutex_index[i].ptr  = mutex;
    }
    onesided_helpmutex_num=0;
    if (MPI_Debug) printf("%5d: In initialize_onesided_helpmutexes: %d mutexes have been initialized.\n",ProcID(),MAX_ONESIDED_HELPMUTEX_ARRAYS);
    }
}

/* Finalization for onesided_helpmutex, free memory for global mutex structure onesided_helpmutex_index */
/* It is called just after finalize_twosided_nxtval */
void finalize_onesided_helpmutexes()
{
    int i;
    if (use_helper_server) {
    /* free the memory for global onesided_helpmutex data structure */
    if (MPI_Debug) printf("%5d: In finalize_onesided_helpmutexes: begin.\n",ProcID());

    for(i=0;i<MAX_ONESIDED_HELPMUTEX_ARRAYS; i++) {
       if ( onesided_helpmutex_index[i].ptr != NULL ) MPIMUTEX_Free(&onesided_helpmutex_index[i].ptr);
    }
    if(onesided_helpmutex_data_struc) {
       if (MPI_Debug) printf("%5d: In finalize_onesided_helpmutexes: free onesided_helpmutex_data_struc.\n",ProcID());
       free(onesided_helpmutex_data_struc);
       onesided_helpmutex_data_struc=NULL;
    }
    if (MPI_Debug) printf("%5d: In finalize_onesided_helpmutexes: %d mutexes have been destroyed.\n",ProcID(),MAX_ONESIDED_HELPMUTEX_ARRAYS);
    }
}

/* allocate one mutex object identified by the original mutex number. Check whether it has overranged, activated, or locked */
int alloc_onesided_helpmutex_orig(int inum)
{
    if (MPI_Debug) printf("%5d: In alloc_onesided_helpmutex_orig begin: original mutex num=%d\n",ProcID(),inum);

    if (inum <0 || inum >= MAX_ONESIDED_HELPMUTEX_ARRAYS ) {
       fprintf(stderr,"ERROR in alloc_onesided_helpmutex_orig: over range! original mutex num=%d, max num=%d\n",inum,MAX_ONESIDED_HELPMUTEX_ARRAYS);
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].actv==1 ) {
       fprintf(stderr,"ERROR in alloc_onesided_helpmutex_orig: original mutex %d has been used by others.\n",inum);
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].lock==1 ) {
       fprintf(stderr,"ERROR in alloc_onesided_helpmutex_orig: original mutex %d has been locked before allocating.\n",inum);
       exit (1);
    }
    else {
       onesided_helpmutex_index[inum].actv=1;
    }

    if (MPI_Debug) printf("%5d: In alloc_onesided_helpmutex_orig end: original mutex %d has been allocated.\n",ProcID(),inum);

    return 0;
}

/* free one mutex object identified by the original mutex number. Check whether it has overranged, activated, or locked */
int free_onesided_helpmutex_orig(int inum)
{

    if (MPI_Debug) printf("%5d: In free_onesided_helpmutex_orig begin: original mutex num=%d\n",ProcID(),inum);

    if (inum <0 || inum >= MAX_ONESIDED_HELPMUTEX_ARRAYS ) {
       fprintf(stderr,"ERROR in free_onesided_helpmutex_orig: over range! original mutex num=%d, max num=%d\n",inum,MAX_ONESIDED_HELPMUTEX_ARRAYS);
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].actv==0 ) {
       fprintf(stderr,"ERROR in free_onesided_helpmutex_orig: original mutex %d has not activated for use.\n",inum);
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].lock==1 ) {
       fprintf(stderr,"WARNING in free_onesided_helpmutex_orig: original mutex %d is still locked before freed. Now unlocking and freeing it...\n",inum);
       unlock_onesided_helpmutex_orig(inum);
       onesided_helpmutex_index[inum].lock = 0;
       onesided_helpmutex_index[inum].actv = 0;
    }
    else {
       onesided_helpmutex_index[inum].actv = 0;
    }

    if (MPI_Debug) printf("%5d: In free_onesided_helpmutex_orig end: original mutex %d has been allocated.\n",ProcID(),inum);

    return 0;
}


/* allocates a set containing the number of mutexes. Only one set of mutexes can exist at a time. Mutexes can be
allocated and freed as many times as needed. Mutexes are numbered: 0, ..., number-1. */
int alloc_onesided_helpmutexes(int number)
{
    int i;
    int handle_orig;

    if (MPI_Debug) printf("%5d: In alloc_onesided_helpmutexes: begin.\n",ProcID());
    if ( number <=0 || number > MAX_ONESIDED_HELPMUTEX_NONGA ) {
       fprintf(stderr,"ERROR in alloc_onesided_helpmutexes: over range! number=%d, max num=%d\n",number,MAX_ONESIDED_HELPMUTEX_NONGA);
       return 1;
    }
    for (i=0;i<number;i++) {
       handle_orig=i+MAX_ONESIDED_HELPMUTEX_GA;
       alloc_onesided_helpmutex_orig(handle_orig);
    }
    onesided_helpmutex_num=number;

    if (MPI_Debug) printf("%5d: In alloc_onesided_helpmutexes end:  %d mutexes have been allocated.\n",ProcID(),number);

    return 0;
}

/* frees the set of mutexes allocated with alloc_onesided_helpmutexes.*/
int free_onesided_helpmutexes()
{
    int i;
    int handle_orig;

    if (MPI_Debug) printf("%5d: In free_onesided_helpmutexes begin: onesided_helpmutex_num=%d\n",ProcID(),onesided_helpmutex_num);

    if (onesided_helpmutex_num<=0) {
       fprintf(stderr,"ERROR in free_onesided_helpmutexes: no mutex is needed to be freed. N mutex=%d\n",onesided_helpmutex_num);
       return 1;
    }
    for(i=0;i<onesided_helpmutex_num; i++) {
       handle_orig=i+MAX_ONESIDED_HELPMUTEX_GA;
       free_onesided_helpmutex_orig(handle_orig);
    }

    onesided_helpmutex_num= 0;

    if (MPI_Debug) printf("%5d: In free_onesided_helpmutexes end: now onesided_helpmutex_num=%d\n",ProcID(),onesided_helpmutex_num);

    return 0;
}

/* lock a mutex object identified by the original mutex number. It is a fatal error for a process
   to attempt to lock a mutex which has already been locked by this process */
int lock_onesided_helpmutex_orig(int inum)
{

    if (MPI_Debug) printf("%5d: In lock_onesided_helpmutex_orig begin: original mutex num=%d\n",ProcID(),inum);

    if (inum <0 || inum >= MAX_ONESIDED_HELPMUTEX_ARRAYS ) {
       fprintf(stderr,"ERROR in lock_onesided_helpmutex_orig: over range! original mutex num=%d, max num=%d\n",inum,MAX_ONESIDED_HELPMUTEX_ARRAYS);
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].actv==0 ) {
       fprintf(stderr,"ERROR in lock_onesided_helpmutex_orig: mutex has not been activated yet.\n");
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].lock==1 ) {
       fprintf(stderr,"ERROR in lock_onesided_helpmutex_orig: attempt to lock a mutex which has already been locked.\n");
       exit (1);
    }
    else {
       MPIMUTEX_Lock(onesided_helpmutex_index[inum].ptr);
       onesided_helpmutex_index[inum].lock = 1;
    }

    if (MPI_Debug) printf("%5d: In lock_onesided_helpmutex_orig: original mutex %d has been locked.\n",ProcID(),inum);

    return 0;
}

/* unlock a mutex object identified by the original mutex number. It is a fatal error for a process
 * to attempt to unlock a mutex which has not been locked by this process. */
int unlock_onesided_helpmutex_orig(int inum)
{

    if (MPI_Debug) printf("%5d: In unlock_onesided_helpmutex_orig begin: original mutex num=%d\n",ProcID(),inum);

    if (inum <0 || inum >= MAX_ONESIDED_HELPMUTEX_ARRAYS ) {
       fprintf(stderr,"ERROR in unlock_onesided_helpmutex_orig: over range! original mutex num=%d, max num=%d\n",inum,MAX_ONESIDED_HELPMUTEX_ARRAYS);
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].actv==0 ) {
       fprintf(stderr,"ERROR in unlock_onesided_helpmutex_orig: mutex has not been activated yet.\n");
       exit (1);
    }
    else if ( onesided_helpmutex_index[inum].lock==0 ) {
       fprintf(stderr,"ERROR in unlock_onesided_helpmutex_orig: attempt to unlock a mutex which has not been locked.\n");
       exit (1);
    }
    else {
       MPIMUTEX_Unlock(onesided_helpmutex_index[inum].ptr);
       onesided_helpmutex_index[inum].lock = 0;
    }

    if (MPI_Debug) printf("%5d: In unlock_onesided_helpmutex_orig: original mutex %d has been unlocked.\n",ProcID(),inum);

    return 0;
}

/* lock a mutex object identified by the wrapped mutex number. It is a fatal error for a process
   to attempt to lock a mutex which has already been locked by this process */
int lock_onesided_helpmutex(int inum)
{
    int inum_orig;

    if (MPI_Debug) printf("%5d: In lock_onesided_helpmutex begin: mutex num=%d, total num=%d\n",ProcID(),inum,onesided_helpmutex_num);

    if (inum <0 || inum >= onesided_helpmutex_num ) {
       fprintf(stderr,"ERROR in mpi_lock_mutex: over range! mutex num=%d, total num=%d\n",inum,onesided_helpmutex_num);
       return 1;
    }
    else {
       inum_orig=inum+MAX_ONESIDED_HELPMUTEX_GA;
       lock_onesided_helpmutex_orig(inum_orig);
    }

    if (MPI_Debug) printf("%5d: In lock_onesided_helpmutex: mutex %d has been locked.\n",ProcID(),inum);

    return 0;
}


/* unlock a mutex object identified by the wrapped mutex number. It is a fatal error for a process
 * to attempt to unlock a mutex which has not been locked by this process. */
int unlock_onesided_helpmutex(int inum)
{
    int inum_orig;

    if (MPI_Debug) printf("%5d: In unlock_onesided_helpmutex begin: mutex num=%d, total num=%d\n",ProcID(),inum,onesided_helpmutex_num);

    if (inum <0 || inum >= onesided_helpmutex_num ) {
       fprintf(stderr,"ERROR in unlock_onesided_helpmutex: over range! mutex num=%d, total num=%d\n",inum,onesided_helpmutex_num);
       return 1;
    }
    else {
       inum_orig=inum+MAX_ONESIDED_HELPMUTEX_GA;
       unlock_onesided_helpmutex_orig(inum_orig);
    }

    if (MPI_Debug) printf("%5d: In unlock_onesided_helpmutex: mutex %d has been unlocked.\n",ProcID(),inum);

    return 0;
}


/* ====================================================================== *\
 * General wrapper routines for both one-sided and two-sided helpmutex    *
 * routines. They are controlled by a switch ( global variable            *
 * use_onesided_helpmutex ).                                              *
\* ====================================================================== */

/* initialization for general helpmutexes -- called just before install_twosided_nxtval */
void initialize_general_helpmutexes()
{
    if ( use_onesided_helpmutex ) initialize_onesided_helpmutexes();
    else initialize_twosided_helpmutexes();
}

/* Finalization for general helpmutex, free memory for global helpmutex structure */
/* It is called just after finalize_twosided_nxtval */
void finalize_general_helpmutexes()
{
    if ( use_onesided_helpmutex ) finalize_onesided_helpmutexes();
    else finalize_twosided_helpmutexes();
}

/* allocate one mutex object identified by the original mutex number. */
int alloc_general_helpmutex_orig(int inum)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = alloc_onesided_helpmutex_orig(inum);
    else mpierr = alloc_twosided_helpmutex_orig(inum);

    return mpierr;
}

/* free one mutex object identified by the original mutex number. */
int free_general_helpmutex_orig(int inum)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = free_onesided_helpmutex_orig(inum);
    else mpierr = free_twosided_helpmutex_orig(inum);

    return mpierr;
}

/* allocates a set containing the number of mutexes. Only one set of mutexes can exist at a time. Mutexes can be
allocated and freed as many times as needed. Mutexes are numbered: 0, ..., number-1. */
int alloc_general_helpmutexes(int number)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = alloc_onesided_helpmutexes(number);
    else mpierr = alloc_twosided_helpmutexes(number);

    return mpierr;
}

/* frees a set of existing mutexes allocated with allocate_general_helpmutexes.*/
int free_general_helpmutexes()
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = free_onesided_helpmutexes();
    else mpierr = free_twosided_helpmutexes();

    return mpierr;
}

/* lock a mutex object identified by the original mutex number. It is a fatal error for a process
   to attempt to lock a mutex which has already been locked by this process */
int lock_general_helpmutex_orig(int inum)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = lock_onesided_helpmutex_orig(inum);
    else mpierr = lock_twosided_helpmutex_orig(inum);

    return mpierr;
}

/* unlock a mutex object identified by the original mutex number. It is a fatal error for a process
 * to attempt to unlock a mutex which has not been locked by this process. */
int unlock_general_helpmutex_orig(int inum)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = unlock_onesided_helpmutex_orig(inum);
    else mpierr = unlock_twosided_helpmutex_orig(inum);

    return mpierr;
}

/* lock a mutex object identified by the wrapped mutex number. It is a fatal error for a process
   to attempt to lock a mutex which has already been locked by this process */
int lock_general_helpmutex(int inum)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = lock_onesided_helpmutex(inum);
    else mpierr = lock_twosided_helpmutex(inum);

    return mpierr;
}

/* unlock a mutex object identified by the wrapped mutex number. It is a fatal error for a process
 * to attempt to unlock a mutex which has not been locked by this process. */
int unlock_general_helpmutex(int inum)
{
    int mpierr;

    if ( use_onesided_helpmutex ) mpierr = unlock_onesided_helpmutex(inum);
    else mpierr = unlock_twosided_helpmutex(inum);

    return mpierr;
}


#else

void mpi_helpmutex_dummy() { }

#endif
