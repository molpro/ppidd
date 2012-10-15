
/* ==================================================================== *\
 * MPI Version of global data structures
\* ==================================================================== */


#ifndef __MPIGA_BASE_H__
#define __MPIGA_BASE_H__

/* ---------------------------------------------------------- *\
   Macro used for global variable -- Ensures unique namespace
\* ---------------------------------------------------------- */
 # define mpigv(a) __mpiga_common__ ## a ## __

/* ------------- *\
***Include Files***
\* ------------- */


 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <strings.h>
 #include <unistd.h>
 #include <time.h>

 #include "mpimutex.h"  /* mpi.h is already included here */
 #include "util/machines.h"

 #include "ppidd_dtype.h"


/* --------------------------------------- *\
   Maximum number of global data structures
   To be safe, make divisible by 8.
\* --------------------------------------- */
 # if !defined MAX_MPI_ARRAYS
 # define MAX_MPI_ARRAYS 512
 # endif

/* MPIGA handle offset */
 # define MPI_GA_OFFSET 1000

/* -------------------------------------------------------- *\
 * Limits on the number of nodes, processors per node, etc. *
\* -------------------------------------------------------- */
 # if !defined MAX_NODES
 # define MAX_NODES 128
 # endif

 # if !defined MAX_SMP_PROCS
 # define MAX_SMP_PROCS 1
 # endif

 # define MAX_PROCESSORS (MAX_NODES*MAX_SMP_PROCS)

 #define MPI2_ONESIDED

/* ------------------------------------------------------------------ *\
   Global Variables -- one-sided array structure and mutex structure
\* ------------------------------------------------------------------ */

#ifdef MPI2_ONESIDED
/* We make MPIGA a pointer to this structure so that users always
   have a pointer, never the actual structure */
typedef struct STRUC_MPIGA {
    char         *name;
    MPI_Win      ga_win;
    void         *win_ptr;
    mpimutex_t   mutex_p;
    /* Datatype and size */
    MPI_Datatype dtype;
    int          dtype_size;
    /* sizes of the global data strcture */
    int          lentot;
    /* number of chunks */
    int          nchunk;
    /* size in each process */
    int          *len;
} *MPIGA;

typedef struct {
/*       short int  irreg;        * 0-regular; 1-irregular distribution  */
/*       int  type;               * data type in array                   */
       int  actv;               /* activity status                      */
       long size;               /* size of data on local process in bytes */
/*       int  elemsize;           * sizeof(datatype)                     */
       MPIGA  ptr;              /* pointer to MPIGA                     */
} mpiglobal_array_t;

   extern mpiglobal_array_t *mpiga_main_data_structure, *MPIGAIndex;

/* mutex structure */
typedef struct {
       int  actv;               /* activity status                      */
       int  lock;               /* lock status (0:unlocked; 1:locked)   */
       mpimutex_t  ptr;         /* pointer to mpimutex_t                */
} mpimutex_t_index;

   /* mpiga_mutex and mpi_helpmutex global data structure */
   extern mpimutex_t_index *mpiga_mutex_data_struc, *mpiga_mutexindex;
   extern mpimutex_t_index *onesided_helpmutex_data_struc, *onesided_helpmutex_index;
#endif


/* ---------------------------------------------------------------------- *\
   Global Variables -- including MPI Global Data Structure and Mutex etc
\* ---------------------------------------------------------------------- */
   extern int mpigv(nga);                      /* Number of MPI Global Data Structure in use */
   extern int mpigv(nmutex);                   /* Number of MPI Mutex in use         */
   extern int onesided_helpmutex_num;          /* Number of one-sided helpmutex in use     */
   extern int use_onesided_helpmutex;          /* Switch of one-sided or two-sided helpmutex in use    */

   extern int mpigv(nprocs);                   /* Number of processes in the group of MPI comm         */
   extern int mpigv(myproc);                   /* Rank of the calling process in the group of MPI comm */

   extern int *mpigv(map);
   /* List of lower and upper indices for local MPIGA that exists on each processor containing a  portion of MPIGA.
    * For a D dimensional MPIGA, the first D elements are the lower indices on 1st processor in proclist,
    * the next D elements are the upper indices on 1st processor in proclist; and then the next D elements are
    * the lower indices on 2nd processor in proclist, and  so on. Here, it very simple since D=1. */
   extern int *mpigv(proclist); /* list of processors containing some portion of MPIGA */

   extern long  mpigv(curmem);    /* currently allocated memory of MPIGA for all processes, exclude freed memory */
   extern long  mpigv(maxmem);    /* maximum limit of memory for MPIGA for all processes */
   extern long  mpigv(grsmem);    /* all allocated gross memory of MPIGA for all processes, include already freed memory */

   extern int NUM_TOTAL_NNODES; /* total nodes for MPI_COMM_WORLD */
   extern int NNODES_SYMMETRY;    /* if all the nodes are symmetric: 1 (yes), 0 (no) */

/* ------------------------------------- *\
   Global Variables -- MPI Communicators
\* ------------------------------------- */
   extern MPI_Comm mpigv(Compute_comm);  /* computation communicator    */



/* =========================== *\
    MPIGA Function Prototypes
\* =========================== */

   int mpiga_initialize(int *, char ***);
   int mpiga_initialize_data(void);
   int mpiga_terminate(void);
   int mpiga_cleanup_finalize(void);
   int mpiga_create_irreg(char *, int *, int , MPI_Datatype , int *);
   int mpiga_create( char *, int , MPI_Datatype , int *);
   int mpiga_free( int );
   long mpiga_localmem(void);
   int mpiga_handle_orig( int );
   int mpiga_inquire_name( int , char **);
   int mpiga_inquire_storetype( int );
   int mpiga_distribution( int , int , int *, int *);
   int mpiga_location( int , int , int , int *, int *, int *);
   int mpiga_put( int , int , int , void *);
   int mpiga_get( int , int , int , void *);
   int mpiga_acc( int , int , int , void * , void *);
   fortint mpiga_read_inc( int ,  int , int );
   int mpiga_zero_patch( int , int , int );
   int mpiga_zero( int );
   int MPI_GSum(MPI_Datatype , void *, int , char *);
   char *strlower(char *);
   int mpiga_create_mutexes(int);
   int mpiga_lock_mutex(int);
   int mpiga_unlock_mutex(int);
   int mpiga_destroy_mutexes(void);

   /* MPI one-sided helpmutex Function Prototypes, from mpi_helpmutex.c */
   void initialize_onesided_helpmutexes(void);
   void finalize_onesided_helpmutexes(void);
   int alloc_onesided_helpmutex_orig(int);
   int free_onesided_helpmutex_orig(int);
   int alloc_onesided_helpmutexes(int);
   int free_onesided_helpmutexes(void);
   int lock_onesided_helpmutex_orig(int);
   int unlock_onesided_helpmutex_orig(int);
   int lock_onesided_helpmutex(int);
   int unlock_onesided_helpmutex(int);

   /* MPI general helpmutex Function Prototypes, from mpi_helpmutex.c */
   void initialize_general_helpmutexes(void);
   void finalize_general_helpmutexes(void);
   int alloc_general_helpmutex_orig(int);
   int free_general_helpmutex_orig(int);
   int alloc_general_helpmutexes(int);
   int free_general_helpmutexes(void);
   int lock_general_helpmutex_orig(int);
   int unlock_general_helpmutex_orig(int);
   int lock_general_helpmutex(int);
   int unlock_general_helpmutex(int);

#endif
