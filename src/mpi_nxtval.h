
/* ==================================================================== *\
 * MPI Version of nextval
 * ============================
 *
\* ==================================================================== */

#ifndef __MPI_NXTVAL_H__
#define __MPI_NXTVAL_H__

/* **Include Files*** */
/*  prerequisite: include mpi.h, and some basic libraries  */
#include <stdint.h>
#include <mpi.h>

/* defined variables */

#define WAKEUPTAG    73333

#define LEN 2
#define INCR 1                     /* increment for NXTVAL */
#define NXTVALFLAG 1555            /* flag of NXTVAL operations: terminate, fetch-and-add, release */
#define COLLECFLAG 1666            /* flag of collective operations for helpga: create, zeroize, destroy */
#define RMAONEFLAG 1777            /* flag of one-element remote memory access operations for integer helpga: get, fetch-and-add, put */
#define RMAETRFLAG 1888            /* flag of extra remote memory access operations for helpga: get/put/accumulate many helpga elements */
#define MUTCOLFLAG 1999            /* flag of collective mutex operations in NXTVAL server: create and destroy mutexes */
#define MUTLOCFLAG 2999            /* flag of lock/unlock mutex operations in NXTVAL server: lock and unlock mutex */

#if !defined MAX_TWOSIDED_HELPGA_ARRAYS
#define MAX_TWOSIDED_HELPGA_ARRAYS 500
#endif
#define TWOSIDED_HELPGA_OFFSET 4000

/* MPI two-sided helpmutex maximum number = (num for GA) + (num for NONGA) */
/* MAX_TWOSIDED_HELPMUTEX_GA=MAX_MPI_ARRAYS, which is defined in mpiga_base.h */
#if !defined MAX_TWOSIDED_HELPMUTEX_GA
  #define MAX_TWOSIDED_HELPMUTEX_GA  500
#endif
#if !defined MAX_TWOSIDED_HELPMUTEX_NONGA
  #define MAX_TWOSIDED_HELPMUTEX_NONGA  100
#endif
#define MAX_TWOSIDED_HELPMUTEX_ARRAYS (MAX_TWOSIDED_HELPMUTEX_GA+MAX_TWOSIDED_HELPMUTEX_NONGA)

/* MPI one-sided helpmutex maximum number = (num for GA) + (num for NONGA) */
#if !defined MAX_ONESIDED_HELPMUTEX_GA
  #define MAX_ONESIDED_HELPMUTEX_GA  100
#endif
#if !defined MAX_ONESIDED_HELPMUTEX_NONGA
  #define MAX_ONESIDED_HELPMUTEX_NONGA  50
#endif
#define MAX_ONESIDED_HELPMUTEX_ARRAYS (MAX_ONESIDED_HELPMUTEX_GA+MAX_ONESIDED_HELPMUTEX_NONGA)
/* MAX_ONESIDED_HELPMUTEX_ARRAYS can't be too large. Otherwise, this error ("too many communicators")
may happen since MPICH2 runs out of context ids. When you create new communicators, it's best to free them after using.
Internally in MPICH2, MPI_Cart_create, MPI_Win_create, and some of the MPI-IO functions call MPI_Comm_dup. */


/* derived data structure */
/* Define a helpga structure STRUC_MPIHELPGA, and pointer to this structure MPIHELPGA.
   Users can always use the pointer to this structure */
typedef struct STRUC_MPIHELPGA {
       char  *name;             /* name of helpga                       */
       void  *ptr_buf;          /* pointer to helpga buffer             */
       MPI_Datatype dtype;      /* Datatype of helpga                   */
       int   nele;              /* total number of elements, but for helper server it is the number of elements located in current server */
       int   *len;              /* size on each compute process         */
       int   *len_help;         /* size on each helper process          */
} *MPIHELPGA;

typedef struct {
       int        actv;              /* activity status                      */
       long       size;              /* size of data "on local process" in bytes */
       MPIHELPGA  ptr;               /* pointer to helpga structure          */
} twosided_helpga_array_t;


typedef struct {
       int        actv;              /* activity status                      */
       int        lock;              /* lock status (0:unlocked; 1:locked)   */
       char       *mutex;            /* pointer to twosided helpmutex buffer, use "char" rather than "int" in order to save space */
} twosided_helpmutex_array_t;


/* ----------------------------------------- *\
   Global Variables -- MPI NXTVAL Data
\* ----------------------------------------- */

    extern twosided_helpga_array_t *twosided_helpga_data_struc, *twosided_helpga_index;
    extern twosided_helpmutex_array_t *twosided_helpmutex_index;            /* twosided helpmutex list */

    extern int *twosided_helpga_map;                  /* List of lower and upper indices for helpga that exists on each helper server containing a portion of helpga */
    extern int *twosided_helpga_proclist;             /* List of serial number of helper server containing some portion of helpga */
    extern int twosided_helpga_num;                   /* Number of MPI helpga in use */
    extern long twosided_helpga_curmem;               /* currently allocated memory of HELPGA for all compute processes, exclude freed memory.
                                                    For helper server it is the memory of HELPGA currently allocated in current server */

    extern int twosided_helpmutex_num;           /* Number of MPI twosided helpmutex in use */
    extern int NPROCS_PER_HELPER;                /* Number of processes(including server) per helper server  */
    extern int NUMBER_OF_SERVER;                 /* total number of servers */


    extern int use_helper_server;                /* helper_server flag: 1 (use); 0 (don't use) */
    extern int SR_parallel;                      /* parallel flag */
    extern MPI_Comm MPIGA_WORK_COMM;             /* MPI Work/Computation Communicator */


/* =============================== *\
    MPI NXTVAL Function Prototypes
\* =============================== */

    extern int TotalNumber_of_Servers();
    extern int NProcs_Work();
    extern int LastServerID();
    extern int SerialNumber_of_Server(int );
    extern int RankNumber_of_Server(int );
    extern int NewRank_of_OldRank(int );
    extern int Server_of_Rank(int );
    extern int Nprocs_of_Server(int );
    extern int OldRank_of_NewRank(int );
    extern void make_worker_comm( MPI_Comm , MPI_Comm * );
    extern void DataHelperServer();
    extern int NXTVAL(int *);
    extern int twosided_helpga_col(int , int );
    extern int twosided_helpga_create_irreg(int , int *, int , int *, char *, int );
    extern int twosided_helpga_create(int , int , int *, char *, int );
    extern int twosided_helpga_locate_server(int , int , int , int *, int *, int *);
    extern int twosided_helpga_distrib( int , int , int *, int *);
    extern int twosided_helpga_location( int , int , int , int *, int *, int *);
    extern int64_t twosided_helpga_one(int , int64_t , int , int);
    extern void twosided_helpga_extra(int , int , int , int , void *);
    extern void twosided_helpga_extra_acc(int , int , int , int , void *, void *);
    extern int twosided_helpga_handle_orig( int );
    extern MPI_Datatype twosided_helpga_inquire_dtype( int );
    extern int twosided_helpga_inquire_name(int , char **);
    extern int twosided_helpmutex_collect(int , int , int );
    extern int twosided_helpmutex_lock(int , int );
    extern int lock_twosided_helpmutex_orig(int );
    extern int unlock_twosided_helpmutex_orig(int );
    extern int lock_twosided_helpmutex(int );
    extern int unlock_twosided_helpmutex(int );
    extern int alloc_twosided_helpmutex_orig(int );
    extern int free_twosided_helpmutex_orig(int );
    extern int alloc_twosided_helpmutexes(int );
    extern int free_twosided_helpmutexes();
    extern void initialize_twosided_helpmutexes();
    extern void finalize_twosided_helpmutexes();
    extern void initialize_twosided_helpga();
    extern int twosided_helpga_release_orig(int );
    extern void finalize_twosided_helpga();
    extern void install_twosided_nxtval();
    extern void finalize_twosided_nxtval();

#endif
/* endif of __MPI_NXTVAL_H__ */
