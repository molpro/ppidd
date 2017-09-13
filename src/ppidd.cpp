#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <climits>
#include <vector>

/*! \file
 * \brief This file contains the PPIDD functions */

#define ga_int int64_t
#define NGA_ACC NGA_Acc64
#define NGA_CREATE NGA_Create64
#define NGA_CREATE_IRREG NGA_Create_irreg64
#define NGA_DISTRIBUTION NGA_Distribution64
#define NGA_LOCATE_REGION NGA_Locate_region64
#define NGA_GET NGA_Get64
#define NGA_PUT NGA_Put64
#define NGA_READ_INC NGA_Read_inc64
#define NGA_ZERO_PATCH NGA_Zero_patch64

#ifdef MPI2
#include "mpiga_base.h"         /* include mpi.h */
#include "mpi_utils.h"
#include "mpi_nxtval.h"
#endif

#ifdef GA_MPI
#include "mpi_utils.h"
#include <ga.h>
#include <ga-mpi.h>
#include <macdecls.h>
#endif

#if defined(MPI2) || defined(GA_MPI)
static int MPIGA_Debug=0;
#endif

#include "ppidd.h"

#ifdef GA_MPI
static int dtype_ga(int dtype) {
 char *errmsg;
 switch (dtype) {
  case 0:
   if (sizeof(FORTINT)==sizeof(int)) return MT_C_INT;
   if (sizeof(FORTINT)==sizeof(long)) return MT_C_LONGINT;
   if (sizeof(FORTINT)==sizeof(long long)) return MT_C_LONGLONG;
   errmsg=strdup(" In dtype_ga: unable to map FORTINT ");
   GA_Error(errmsg,dtype);
   free(errmsg);
   break;
  case 1:
   return MT_C_DBL;
  default:
   errmsg=strdup(" In dtype_ga: wrong data type ");
   GA_Error(errmsg,dtype);
   free(errmsg);
 }
 return -1;
}

static size_t dtype_ga_bytes(int dtype) {
 switch (dtype) {
  case 0:
   return sizeof(FORTINT);
  case 1:
   return sizeof(double);
  default:
   char *errmsg=strdup(" In dtype_ga_bytes: wrong data type ");
   GA_Error(errmsg,dtype);
   free(errmsg);
 }
 return -1;
}
#endif

extern "C" {

/*! Initialize the PPIDD parallel environment
 *
 *  - For \b GA, includes initialization of MPI and GA.
 *  - For \b MPI2, calls MPI_Init.
 */
   void PPIDD_Initialize(int argc, char **argv) {
#ifdef MPI2
    int mpierr=mpiga_initialize(&argc,&argv);
    mpi_test_status("PPIDD_Initialize:",mpierr);
#elif defined(GA_MPI)
    MPI_Init(&argc, &argv);                     /* initialize MPI */
    GA_Initialize_args(&argc,&argv);            /* initialize GA */
#endif
   }


/*! Initialize the PPIDD data structure
 *
 *  - For \b GA, does nothing
 *  - For \b MPI2, Initialize global data structure and set helper server.
 */
   void PPIDD_Initialize_data(void) {
#ifdef MPI2
      mpiga_initialize_data();
      if(MPIGA_Debug)printf("%5d: [PPIDD_Initialize_data] end.\n",ProcID());
#endif
   }


/*! Return communicator for worker process group, which excludes any helper process.
 * For simplicity, the returned data type is set as integer. So for C calling, it may have to be converted
 * into MPI_Comm by MPI_Comm_f2c afterward if needed.
 *
 *  - For \b MPI2, Return communicator for worker process group.
 *  - For \b GA and serial cases, should not be called.
 */
   int64_t PPIDD_Worker_comm(void) {
#ifdef MPI2
      MPI_Comm mycomm=mpigv(Compute_comm);
      MPI_Fint fcomm;

/* test whether worker communicator contains all the processes, if so then return MPI_COMM_WORLD */
      int np_all, np_worker=mpigv(nprocs);
      MPI_Comm_size(MPI_COMM_WORLD, &np_all);
      if(np_all==np_worker) mycomm=MPI_COMM_WORLD;

      fcomm=MPI_Comm_c2f(mycomm);
      return (int64_t)fcomm;
#elif defined(GA_MPI)
      MPI_Comm mpicomm = GA_MPI_Comm();
      MPI_Fint fcomm=MPI_Comm_c2f(mpicomm);
      return (int64_t)fcomm;
#else
      fprintf(stderr," ERROR: PPIDD_Worker_comm should not be called in NON-MPI2 cases.\n");
      return (int64_t)-1;
#endif
   }


/*! Terminate the PPIDD parallel environment.
 *
 *  - For \b GA, tidy up global arrays and MPI, analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#TERMINATE
 *  - For \b MPI2, tidy up some associated resources and call MPI_Finalize.
 */
   void PPIDD_Finalize(void) {
#ifdef MPI2
      mpiga_terminate();
#elif defined(GA_MPI)
      GA_Terminate();
      MPI_Finalize();
#endif
   }


/*! \cond */
/*  Detect whether MA is used for allocation of GA memory.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#USES_MA
 *  - \b MPI2 always returns <tt>.false.</tt>
 */
   int PPIDD_Uses_ma() {
#ifdef MPI2
      return 0;
#elif defined(GA_MPI)
      if (GA_Uses_ma()) return 1;
      else return 0;
#else
      return 0;
#endif
  }


/*  Initialize the memory allocator.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/ma/MAapi.html
 *  - \b MPI2 always returns <tt>.true.</tt>
 */
   int PPIDD_MA_init(int dtype, int64_t *stack, int64_t *heap) {
#ifdef GA_MPI
      Integer istack=(Integer)*stack;
      Integer iheap=(Integer)*heap;
      Integer gadtype=(Integer)dtype_ga(dtype);
      if( MA_init(gadtype, istack, iheap)) return 1;
      else return 0;
#else
      return 1;
#endif
  }
/*! \endcond */

/*! Return an elapsed time on the calling process.
 *
 * Returns a floating-point number of seconds, representing elapsed wall-clock time since an arbitrary time in the past.
 * The times returned are local to the node that called them. There is no requirement that different nodes return the same time.
 * This is a local operation.
 *
 *  - \b GA calls GA_Wtime, which is available only in release 4.1 or greater, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#WTIME
 *  - \b MPI2 calls MPI_Wtime
 */
   void PPIDD_Wtime(double *ctime) {
#ifdef MPI2
      *ctime = MPI_Wtime();
#elif defined(GA_MPI)
      *ctime = GA_Wtime();
#else
      *ctime = (double)0;
#endif
  }


/*! Print an error message and abort the program execution.
 *
 *  - \b GA calls GA_Error, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ERROR
 *  - For \b MPI2, prints error, and then calls MPI_Abort.
 */
   void PPIDD_Error(char *message,int *code) {
#ifdef MPI2
      MPIGA_Error(message,*code);
#elif defined(GA_MPI)
      GA_Error(message,*code);
#else
      fprintf(stdout," %s %d (%#x).\n", message,*code,*code);
      fflush(stdout);
      fprintf(stderr," %s %d (%#x).\n", message,*code,*code);

      printf(" PPIDD_Error: now exiting...\n");
      exit(1);
#endif
   }


/*! Set the flag for data helper server and set the number of how many processes own one helper server
 *
 *  Set the helper_server flag: 1 (use); 0 (don't use). */
   void PPIDD_Helper_server(int *flag, int64_t *numprocs_per_server) {
#ifdef MPI2
      if ((int)*flag) {                              /* mutilple helper servers, node helper server, and single helper server */
         use_helper_server=1;
         if ( (int)*numprocs_per_server > 1 ) {      /* reasonable mutilple helper servers, and single helper server */
            NPROCS_PER_HELPER=(int)*numprocs_per_server;
         }
         else if ( (int)*numprocs_per_server == 0 ) {/* node helper server: one helper server on every node */
           int mpinp;
           MPI_Comm_size(MPI_COMM_WORLD, &mpinp);
           if (NNODES_SYMMETRY) {
             NPROCS_PER_HELPER=mpinp/NUM_TOTAL_NNODES;
             if (NPROCS_PER_HELPER==1) {             /* node helper server: if NPROCS_PER_HELPER==1, then use only one single helper server */
               if(NUM_TOTAL_NNODES>1)fprintf(stdout,"%5d: WARNING: only one process on each node. Will use only one single helper server for all processes!\n", ProcID());
               NPROCS_PER_HELPER=99999999;
             }
           }
           else {                                    /* node helper server: if not all the nodes are symmetric, then use only one single helper server */
             fprintf(stderr,"%5d: WARNING: not all the nodes are symmetric. Will use only one single helper server for all processes!\n", ProcID());
             NPROCS_PER_HELPER=99999999;
           }
         }
         else {                                      /* unreasonable mutilple helper servers, then only use single helper server */
            fprintf(stderr,"%5d: WARNING: nprocs_per_server=%d is unreasonable. Will use only single helper server for all processes!\n", ProcID(),(int)*numprocs_per_server);
            NPROCS_PER_HELPER=99999999;
         }
      }
      else {                                         /* no helper server */
         use_helper_server=0;
         NPROCS_PER_HELPER=-1;
      }
#endif
   }


/*!  Determine the total number of processes available (including helper process if there is one).
 *
 *  - \b GA calls GA_Nnodes,
 *  http://hpc.pnl.gov/globalarrays/api/c_op_api.html#NNODES
 *  - \b MPI2 calls MPI_Comm_size for communicator MPI_COMM_WORLD.
 */
   void PPIDD_Size_all(int64_t *np) {
#ifdef MPI2
      int mpinp;
      MPI_Comm mpicomm=MPI_COMM_WORLD;

      MPI_Comm_size(mpicomm, &mpinp);
      *np = (int64_t) mpinp;
#elif defined(GA_MPI)
      *np = (int64_t)GA_Nnodes();
#else
      *np = (int64_t)1;
#endif
   }


/*  Determine the number of compute processes.
 *
 *  - \b GA calls GA_Nnodes,
 *  http://hpc.pnl.gov/globalarrays/api/c_op_api.html#NNODES
 *  - \b MPI2 calls MPI_Comm_size for computational communicator.
 */
   void PPIDD_Size(int64_t *np) {
#ifdef MPI2
      int mpinp;
      MPI_Comm mpicomm=mpigv(Compute_comm);

      MPI_Comm_size(mpicomm, &mpinp);
      *np = (int64_t) mpinp;
#elif defined(GA_MPI)
      *np = (int64_t)GA_Nnodes();
#else
      *np = (int64_t)1;
#endif
   }


/*! Determine the rank of the calling process.
 *
 *  - \b GA calls GA_Nodeid,
 *  http://hpc.pnl.gov/globalarrays/api/c_op_api.html#NODEID
 *  - \b MPI2 calls MPI_Comm_rank in computational communicator.
 */
   void PPIDD_Rank(int64_t *me) {
#ifdef MPI2
      int mpime;
      MPI_Comm mpicomm=mpigv(Compute_comm);

      MPI_Comm_rank(mpicomm, &mpime);
      *me = (int64_t) mpime;
#elif defined(GA_MPI)
      *me = (int64_t)GA_Nodeid();
#else
      *me = (int64_t)0;
#endif
   }


/*! Initialize tracing of completion status of data movement operations.
 *
 *  - \b GA calls GA_Init_fence, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#INIT_FENCE
 *  - \b MPI2 does nothing
 */
   void PPIDD_Init_fence(void) {
#ifdef GA_MPI
      GA_Init_fence();
#endif
   }


/*! Block the calling process until all data transfers complete.
 *
 *  - \b GA calls GA_Fence, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#FENCE
 *  - \b MPI2 does nothing
 */
   void PPIDD_Fence(void) {
#ifdef GA_MPI
      GA_Fence();
#endif
   }



#if defined(MPI2) || defined(GA_MPI)
/* =================== nonblocking mpi message list ================= */
#define MAX_MPIQ_LEN 128      /* Maximum no. of outstanding messages */
/* In application programs, nonblocking send/recv may be never used, so here set it to a small number. It could be increased if necessary. */
/*! \cond */
static struct msg_mpiq_struct{
  MPI_Request request;
  long         node;
  long         type;
  long         lenbuf;
  long         snd;
  long         from;
} msg_mpiq[MAX_MPIQ_LEN];
/*! \endcond */

static int n_in_msg_mpiq=0;
/* =================================================================== */
#endif


/*! Blocking/nonblocking send.
 *
 *  - \b GA calls SND_.
 *  - \b MPI2 calls MPI_Send ( sync is 1) or MPI_Isend ( sync is 0).
 */
   void PPIDD_Send(void *buf,int64_t *count,int dtype,int64_t *dest,int64_t *sync) {
#if defined(MPI2) || defined(GA_MPI)
  #ifdef MPI2
      MPI_Comm mpicomm=mpigv(Compute_comm);
  #endif
  #ifdef GA_MPI
      MPI_Comm mpicomm = GA_MPI_Comm();
  #endif
      int mpicount=(int)*count;
      int mpidest=(int)*dest;
      int mpitag=dtype;
      int mpisync=(int)*sync;
      int mpilenbuf,mpierr;
      MPI_Datatype mpidtype;
      int sizempidtype;

      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      mpilenbuf=mpicount*sizempidtype;

      if (MPIGA_Debug) {
         printf("PPIDD_SEND: node %d sending to %d, len(bytes)=%d, mes tag=%d, sync=%d\n",
                 ProcID(), mpidest, mpilenbuf, mpitag, mpisync);
         fflush(stdout);
      }

      if (mpisync) {
         mpierr=MPI_Send(buf,mpilenbuf,MPI_CHAR,mpidest,mpitag,mpicomm);
         mpi_test_status("PPIDD_SEND: SEND:",mpierr);
      }
      else {
         if (n_in_msg_mpiq >= MAX_MPIQ_LEN) {
            MPIGA_Error("PPIDD_SEND: nonblocking SEND: overflowing async Queue limit",n_in_msg_mpiq);
         }
         mpierr = MPI_Isend(buf, mpilenbuf, MPI_CHAR,mpidest, mpitag,mpicomm,
                     &msg_mpiq[n_in_msg_mpiq].request);
         mpi_test_status("PPIDD_SEND: nonblocking SEND:",mpierr);

         msg_mpiq[n_in_msg_mpiq].node   =(long) mpidest;
         msg_mpiq[n_in_msg_mpiq].type   =(long) dtype;
         msg_mpiq[n_in_msg_mpiq].lenbuf =(long) mpilenbuf;
         msg_mpiq[n_in_msg_mpiq].snd = (long)1;
      }
#else
      printf(" ERROR: PPIDD_Send should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Blocking/nonblocking receive.
 *
 *  - \b GA calls RCV_.
 *  - \b MPI2 calls MPI_Recv ( sync is 1) or MPI_Irecv( sync is 0).
 */
   void PPIDD_Recv(void *buf,int64_t *count,int dtype,int64_t *source,int64_t *lenreal,int64_t *sourcereal,int64_t *sync) {
#if defined(MPI2) || defined(GA_MPI)
  #ifdef MPI2
      MPI_Comm mpicomm=mpigv(Compute_comm);
  #endif
  #ifdef GA_MPI
      MPI_Comm mpicomm = GA_MPI_Comm();
  #endif
      int mpicount=(int)*count;
      int mpitag=dtype;
      int mpisource=(int)*source;
      int mpisync=(int)*sync;
      int mpinode,mpilenbuf,mpierr;
      MPI_Status status;
      MPI_Request request;
      MPI_Datatype mpidtype;
      int sizempidtype;

      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      mpilenbuf=mpicount*sizempidtype;

      if (mpisource == -1)
         mpinode = MPI_ANY_SOURCE;
      else
         mpinode = mpisource;

      if (MPIGA_Debug) {
         printf("PPIDD_Recv: node %d receving from %d, len(bytes)=%d, mes tag=%d, sync=%d\n",
                 ProcID(), mpisource, mpilenbuf, mpitag, mpisync);
         fflush(stdout);
      }

      if(mpisync==0){
         if (n_in_msg_mpiq >= MAX_MPIQ_LEN) {
            MPIGA_Error("PPIDD_Recv: nonblocking RECV: overflowing async Queue limit", n_in_msg_mpiq);
         }
         mpierr = MPI_Irecv(buf,mpilenbuf,MPI_CHAR,mpinode,mpitag,mpicomm,&request);
         mpi_test_status("PPIDD_Recv: nonblocking RECV:",mpierr);

         *sourcereal = (int64_t) mpinode;          /* Get source node  */
         *lenreal = (int64_t) (-1);
         msg_mpiq[n_in_msg_mpiq].request = request;
         msg_mpiq[n_in_msg_mpiq].node   = (long)*source;
         msg_mpiq[n_in_msg_mpiq].type   = (long)dtype;
         msg_mpiq[n_in_msg_mpiq].lenbuf = (long)mpilenbuf;
         msg_mpiq[n_in_msg_mpiq].snd = (long)0;
         n_in_msg_mpiq++;
      }
      else{
         mpierr = MPI_Recv(buf,mpilenbuf,MPI_CHAR,mpinode,mpitag,mpicomm,&status);
         mpi_test_status("PPIDD_RECV: RECV:",mpierr);
         mpierr = MPI_Get_count(&status, MPI_CHAR, &mpilenbuf);
         mpi_test_status("PPIDD_RECV: Get_count:",mpierr);
         *sourcereal = (int64_t)status.MPI_SOURCE;
         *lenreal    = (int64_t)mpilenbuf;
      }
#else
      printf(" ERROR: PPIDD_Recv should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Wait for completion of all asynchronous send/receive.
 *
 * ignores nodesel !! */
/*!
 *  - \b MPI2/GA_MPI calls MPI_Wait for all asynchronous requests.
 */
   void PPIDD_Wait(int64_t *nodesel) {
#if defined(MPI2) || defined(GA_MPI)
      int mpierr,i;
      MPI_Status status;

      for (i=0; i<n_in_msg_mpiq; i++){
         if (MPIGA_Debug) {
            printf("PPIDD_Wait: node %d waiting for msg to/from %ld, #%d\n", ProcID(), msg_mpiq[i].node, i);
            fflush(stdout);
         }
         mpierr = MPI_Wait(&msg_mpiq[i].request, &status);
         mpi_test_status("PPIDD_Wait:",mpierr);
      }
      n_in_msg_mpiq = 0;
#endif
   }


/*! Detect whether a message of the specified type is available.
 *
 *  Return <tt>.true.</tt> if the message is available, otherwise <tt>.false.</tt>.
 *
 *  - \b MPI2/GA_MPI calls MPI_Iprobe
 */
   int PPIDD_Iprobe(int64_t *tag,int64_t *source) {
#if defined(MPI2) || defined(GA_MPI)
  #ifdef MPI2
      MPI_Comm mpicomm=mpigv(Compute_comm);
  #endif
  #ifdef GA_MPI
      MPI_Comm mpicomm = GA_MPI_Comm();
  #endif
      int mpitag=(int)*tag;
      int mpisource,mpierr,flag;
      MPI_Status status;

      mpisource = (*source < 0) ? MPI_ANY_SOURCE  : (int) *source;
      mpierr = MPI_Iprobe(mpisource, mpitag, mpicomm, &flag, &status);
      mpi_test_status("PPIDD_Iprobe:",mpierr);
      if(flag) return 1 ;
      else return 0 ;
#else
    return 0 ;
#endif
   }


/*! Broadcast a message from the root process to all other processes.
 *
 *  Collective operation.
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#BRDCST
 *  - \b MPI2 calls MPI_Bcast
 *
 *  - \c dtype=0 : Fortran integer and logical types
 *  - \c dtype=1 : Fortran double precision type */
   void PPIDD_BCast(void *buffer,int64_t *count,int dtype,int64_t *root) {
#ifdef MPI2
      int mpicount=(int)*count;
      int mpiroot=(int)*root;
      MPI_Datatype mpidtype;
      int sizempidtype;
      MPI_Comm mpicomm;
      int mpierr;

      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      mpicomm=mpigv(Compute_comm);

      mpierr=MPI_Bcast(buffer,mpicount,mpidtype,mpiroot,mpicomm);
      mpi_test_status("PPIDD_BCast:",mpierr);
#elif defined(GA_MPI)
      int gacount=(int)*count;
      int garoot=(int)*root;
      size_t ctype=dtype_ga_bytes(dtype);
      int galenbuf=ctype*gacount;
      GA_Brdcst(buffer, galenbuf, garoot);
#endif
   }


/*! Synchronize processes and ensure all have reached this routine.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#SYNC
 *  - \b MPI2 calls MPI_Barrier
 */
   void PPIDD_Barrier(void) {
#ifdef MPI2
      int mpierr;

      mpierr=MPI_Barrier(mpigv(Compute_comm));
      mpi_test_status("PPIDD_Barrier:",mpierr);
#elif defined(GA_MPI)
      GA_Sync();
#endif
   }



/*! Combine values from all processes and distribute the result back to all processes.
 *
 *  - \b GA analogous to GA_dgop and GA_igop, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#GOP
 *  - For \b MPI2, it is realized by MPI_Allreduce
 *
 *  - \c type=0 : Fortran Integer
 *  - \c type=1 : Fortran Double Precision */
   void PPIDD_Gsum(int dtype,void *buffer,int64_t *len, char *op) {
#ifdef MPI2
      MPI_Datatype mpidtype;
      int sizempidtype;
      int mpilen=(int)*len;
      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      MPI_GSum(mpidtype,buffer,mpilen, op);
#elif defined(GA_MPI)
      int buflen=(int)*len;
      int gadtype=dtype_ga(dtype);
      GA_Gop(gadtype, buffer, buflen, op);
#endif
   }




/*! Create an array by following the user-specified distribution and return integer handle representing the array.
 *
 * Irregular distributed array data are stored across the distributed processes.
 *  - \c dtype=0 : Fortran integer and logical types
 *  - \c dtype=1 : Fortran double precision type

 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#CREATE_IRREG
 */
   int PPIDD_Create_irreg(char *name, int64_t *lenin, int64_t *nchunk, int dtype, int64_t *storetype, int64_t *handle) {
#ifdef MPI2
      int mpierr;
      int mpinchunk=(int)*nchunk;
      int stype=(int)*storetype;
      MPI_Datatype mpidtype;
      int sizempidtype;
      int mpihandle;

      std::vector<int> mpilenin(mpinchunk);
      for (int i=0;i<mpinchunk;i++) mpilenin[i]=(int)lenin[i];
      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      if (use_helper_server==0) {
        mpierr=mpiga_create_irreg(name, &mpilenin[0], mpinchunk, mpidtype, &mpihandle);
      }
      else {
        if (stype==0)
          mpierr=mpiga_create_irreg(name, &mpilenin[0], mpinchunk, mpidtype, &mpihandle);
        else {
          int mproc=0;
          mpierr=twosided_helpga_create_irreg(mproc, &mpilenin[0], mpinchunk, &mpihandle, name, mpidtype);
        }
      }
      *handle=(int64_t)mpihandle;
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int ndim=1;
      ga_int nblock=(ga_int)*nchunk;
      int np;
      int i;
      ga_int iad;
      int gahandle;
      int gadtype=dtype_ga(dtype);

      ga_int block[1]={nblock};
      np = GA_Nnodes();
/* map[np] or map[nblock] ? */
      std::vector<ga_int> map(np);

      for(iad=0,i=0;i<nblock;i++){
        map[i]=iad;
        iad=iad+(ga_int)lenin[i];
      }
      for(i=nblock;i<np;i++) map[i]=iad;
      ga_int dims[1]={iad};

/*      printf("\n NGA_CREATE_IRREG: %s created, dims=%d, ndim=%d\n",name,dims[1],ndim); */
      gahandle=NGA_CREATE_IRREG(gadtype, ndim, dims, name, block, &map[0]);

      *handle=(int64_t)gahandle;
      return 1 ;
#else
      printf(" ERROR: PPIDD_Create_irreg should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Create an array using the regular distribution model and return integer handle representing the array.
 *
 *  - \c dtype=0  : Fortran integer and logical types
 *  - \c dtype=1  : Fortran double precision type
 *  - \c storetype=0 : Normal distributed array stored across the distributed processes
 *  - \c storetype>=1: Low-latency array stored on one or more helper processes (effective only when helper process is enabled). \c storetype is advisory: the underlying implementation will use up to \c storetype helpers.

 *  - For \b GA, storetype doesn't take effect, and data are always stored across the distributed processes, analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#CREATE
 *  - For \b MPI2, the library can presently be built with zero or one (default) helpers.
 *       When helper process is disabled, \c storetype doesn't take effect, and data are always stored across the distributed processes.
 */
   int PPIDD_Create(char *name,int64_t *lentot, int dtype, int64_t *storetype, int64_t *handle) {
#ifdef MPI2
      int mpierr;
      int mpilentot;
      MPI_Datatype mpidtype;
      int sizempidtype;
      int stype=(int)*storetype;
      int mpihandle;

      if (*lentot > INT_MAX) {
       printf(" ERROR: PPIDD_Create: lentot too large for MPI\n");
       exit(1);
      }
      else mpilentot=(int)*lentot;
      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      if (use_helper_server==0) {
        mpierr=mpiga_create( name, mpilentot, mpidtype, &mpihandle );
      }
      else {
        if (stype==0)
         mpierr=mpiga_create( name, mpilentot, mpidtype, &mpihandle );
        else {
         int mproc=0;
         mpierr=twosided_helpga_create(mproc, mpilentot, &mpihandle, name, mpidtype);
        }
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Create: array %s created, dtype=%d, storetype=%d\n",ProcID(),name,stype,dtype);

      *handle=(int64_t)mpihandle;
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      ga_int galentot=(ga_int)*lentot;
      int gahandle;
      int gadtype=dtype_ga(dtype);

      ga_int dims[1]={galentot};
      ga_int block[1]={-1};

/*      printf("\n NGA_CREATE: %s created, dims=%d, ndim=%d\n",name,*dims,ndim); */
      gahandle=NGA_CREATE(gadtype, 1, dims, name, block);

      *handle=(int64_t)gahandle;
      return 1 ;
#else
      printf(" ERROR: PPIDD_Create should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Deallocate the array represented by handle and free any associated resources.
 *
 *  - \b GA analogous http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DESTROY
 */
   int PPIDD_Destroy(int64_t *handle) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_free(mpihandle);
      else {
         int mproc=-NProcs_Work();
         mpierr=twosided_helpga_col(mproc, mpihandle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Destroy: array %d destroyed!\n",ProcID(),mpihandle);
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int ihandle = (int) *handle;
      GA_Destroy(ihandle);
      return 1 ;
#else
      printf(" ERROR: PPIDD_Destroy should not be called in serial case.\n");
      exit(1);
#endif
   }



/*! Return the range of a distributed array held by a specified process.
 *
 *  Return <tt>.true.</tt> if successful, otherwise <tt>.false.</tt>
 *  If no array elements are owned by the process, the range is returned as [0,-1].
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DISTRIBUTION
 */
   int PPIDD_Distrib(int64_t *handle,int64_t *rank,int64_t *ilo,int64_t *ihi) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpirank=(int)*rank;
      int mpiilo;
      int mpiihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 ) {
         mpierr=mpiga_distribution( mpihandle, mpirank, &mpiilo, &mpiihi);
      }
      else {
         mpierr=twosided_helpga_distrib( mpihandle, mpirank, &mpiilo, &mpiihi);
      }

      *ilo = (int64_t) mpiilo;
      *ihi = (int64_t) mpiihi;
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      int garank=(int)*rank;
      ga_int gailo[1];
      ga_int gaihi[1];

      NGA_DISTRIBUTION(gahandle, garank, gailo, gaihi);
/* If no array elements are owned by process iproc, the range is returned as lo[ ]=0 and hi[ ]= -1 for all dimensions. */
      if (gailo[0]<=gaihi[0]) {
         *ilo = (int64_t) (gailo[0] + 1);
         *ihi = (int64_t) (gaihi[0] + 1);
      }
      else {
         *ilo = (int64_t) (gailo[0]);
         *ihi = (int64_t) (gaihi[0]);
      }
      return 1 ;
#else
      printf(" ERROR: PPIDD_Distrib should not be called in serial case.\n");
      exit(1);
#endif
   }



/*! Return a list of the processes that hold the data.
 *
 *  Parts of the specified patch might be actually 'owned' by several processes.
 *  np is the number of processes hold tha data (return 0  if ilo/ihi are out of bounds "0").
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#LOCATE_REGION
 */
   int PPIDD_Location(int64_t *handle,     /*!< array handle */
                       int64_t *ilo,        /*!< lower element subscript, 1 (not 0) for fisrt element */
                       int64_t *ihi,        /*!< higher element subscript */
                       int64_t *map,        /*!< return start/end index for \c proclist */
                       int64_t *proclist,   /*!< proc id list */
                       int64_t *np          /*!< proc number */
   ) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpisize,mpinp;
      int mpierr;
      MPI_Comm mpicomm;

      mpicomm=mpigv(Compute_comm);
      MPI_Comm_size(mpicomm, &mpisize);
      std::vector<int> mpimap(2*mpisize);
      std::vector<int> mpiproclist(mpisize);

      if ( mpiga_inquire_storetype(mpihandle) == 0 ) {
         mpierr=mpiga_location( mpihandle, mpiilo, mpiihi, &mpimap[0], &mpiproclist[0], &mpinp);
      }
      else {
         mpierr=twosided_helpga_location( mpihandle, mpiilo, mpiihi, &mpimap[0], &mpiproclist[0], &mpinp);
      }

      for (int i=0;i<mpinp;i++) {
         map[2*i]=(int64_t)mpimap[2*i];
         map[2*i+1]=(int64_t)mpimap[2*i+1];
         proclist[i]=(int64_t)mpiproclist[i];
      }
      *np = (int64_t) mpinp;
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int mpihandle=(int)*handle;
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      int mpisize = GA_Nnodes();
      std::vector<ga_int> mpimap(2*mpisize);
      std::vector<int> mpiproclist(mpisize);

      int mpinp=NGA_LOCATE_REGION( mpihandle, mpiilo, mpiihi, &mpimap[0], &mpiproclist[0]);

      for (int i=0;i<mpinp;i++) {
         map[2*i]=(int64_t)(mpimap[2*i]+1);
         map[2*i+1]=(int64_t)(mpimap[2*i+1]+1);
         proclist[i]=(int64_t)mpiproclist[i];
      }
      *np = (int64_t) mpinp;
      return 1 ;
#else
      printf(" ERROR: PPIDD_Location should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Copies data from array section to the local array buffer according to starting and ending index.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#GET
 */
   int PPIDD_Get(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_get(mpihandle, mpiilo, mpiihi, buff);
      else {
         int mproc=0;
         int ielem=mpiilo;
         int64_t val;
         MPI_Datatype dtype=twosided_helpga_inquire_dtype(mpihandle);
         if ( (mpiilo==mpiihi) && (dtype==MPI_INT||dtype==MPI_LONG||dtype==MPI_LONG_LONG) ) { /* PPIDD_helpga_get_inum */
            int64_t nelem_valput=1;
            int64_t *ibuff;

            val=twosided_helpga_one(mproc, nelem_valput, ielem, &mpihandle);
            ibuff=(int64_t *)buff;
            *ibuff=(int64_t)val;
         }
         else if (mpiilo <= mpiihi) { /* PPIDD_helpga_get */
            int nelem=mpiihi-mpiilo+1;

            val=twosided_helpga_extra(mproc, nelem, ielem, &mpihandle,buff);
         }
         else {
            MPIGA_Error("PPIDD_Get: starting index > ending index, handle=",mpihandle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Get: Get value from array handle= %d [%d--%d].\n",ProcID(),mpihandle,mpiilo,mpiihi);
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int mpihandle=(int)*handle;
      ga_int ld[1]={1};
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_GET(mpihandle, mpiilo, mpiihi, buff, ld);
      return 1 ;
#else
      printf(" ERROR: PPIDD_Get should not be called in serial case.\n");
      exit(1);
#endif
   }

/*! Put local buffer data into a section of a global array according to starting and ending index.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#PUT
 */
   int PPIDD_Put(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_put(mpihandle, mpiilo, mpiihi, buff);
      else {
         int mproc=-NProcs_Work();
         int ielem=mpiilo;
         MPI_Datatype dtype=twosided_helpga_inquire_dtype(mpihandle);
         if ( (mpiilo==mpiihi) && (dtype==MPI_INT||dtype==MPI_LONG||dtype==MPI_LONG_LONG) ) { /* PPIDD_helpga_put_inum */
            int64_t *ibuff=(int64_t *)buff;
            int64_t valput=(int64_t)*ibuff;
            twosided_helpga_one(mproc, valput, ielem, &mpihandle);
         }
         else if (mpiilo <= mpiihi) { /* PPIDD_helpga_put */
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra(mproc, nelem, ielem, &mpihandle,buff);
         }
         else {
            MPIGA_Error("PPIDD_Put: starting index > ending index, handle=",mpihandle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Put: Put buff numbers to array handle=%d [%d--%d].\n",ProcID(),mpihandle,mpiilo,mpiihi);
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int mpihandle=(int)*handle;
      ga_int ld[1]={1};
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_PUT(mpihandle, mpiilo, mpiihi, buff, ld);
      return 1 ;
#else
      printf(" ERROR: PPIDD_Put should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Accumulate data into a section of a global array.
 *
 * Atomic operation.  global array section (ilo, ihi) += *fac * buffer
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ACC
 */
   int PPIDD_Acc(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_acc(mpihandle, mpiilo, mpiihi, buff, fac);
      else {
         if (mpiilo <= mpiihi) { /* PPIDD_helpga_acc */
            int mproc=NProcs_Work();
            int ielem=mpiilo;
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra_acc(mproc, nelem, ielem, &mpihandle, buff, fac);
         }
         else {
            MPIGA_Error("PPIDD_Put: starting index > ending index, handle=",mpihandle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Acc: Accumulate buff numbers to array handle=%d [%d--%d].\n",ProcID(),mpihandle,mpiilo,mpiihi);
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int mpihandle=(int)*handle;
      ga_int ld[1]={1};
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_ACC(mpihandle, mpiilo, mpiihi, buff, ld, fac);
      return 1 ;
#else
      printf(" ERROR: PPIDD_Acc should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Atomically read and increment an element in an integer array.
 *
 *  Reads data from the (inum) element of a global array of integers, returns that value, and increments the (inum)
    element by a given increment. This is a fetch-and-add operation.
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#READ_INC
 */
   void PPIDD_Read_inc(int64_t *handle,int64_t *inum,int64_t *incr,int64_t *returnval) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      int mpiinum = (int) *inum;
      int mpiincr = (int) *incr;
      int64_t mpivalue;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpivalue=(int64_t)mpiga_read_inc(mpihandle,mpiinum,mpiincr);
      else {                                              /* PPIDD_helpga_readinc */
         int mproc=NProcs_Work();

         int64_t vincr=(int64_t)*incr;
         mpivalue=twosided_helpga_one(mproc, vincr, mpiinum, &mpihandle);
      }
      *returnval=(int64_t)mpivalue;
      if(MPIGA_Debug)printf("%5d: In PPIDD_Read_inc: fetch-and-add element[%d] of array handle=%d by increment=%d\n",
                            ProcID(),mpiinum,mpihandle,mpiincr);
#elif defined(GA_MPI)
      int ihandle = (int) *handle;
      ga_int mpiinum[1];
      long gaincr = (long) *incr;
      long gavalue;

      mpiinum[0] = (ga_int) *inum-1;
      gavalue=NGA_READ_INC(ihandle,mpiinum, gaincr);
      *returnval=(int64_t)gavalue;
#else
      printf(" ERROR: PPIDD_Read_inc should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Set all the elements in an array patch to zero.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ZERO_PATCH
 */
   void PPIDD_Zero_patch(int64_t *handle,int64_t *ilo,int64_t *ihi) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      int mpiilo = (int) *ilo;
      int mpiihi = (int) *ihi;
      int mpierr=0;
      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_zero_patch(mpihandle,mpiilo,mpiihi);
      else
         MPIGA_Error("PPIDD_Zero_patch: invalid storetype, should be 0. handle=",mpihandle);

      if(mpierr!=0) MPI_Abort(mpigv(Compute_comm),911);
#elif defined(GA_MPI)
      int ihandle = (int) *handle;
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_ZERO_PATCH(ihandle, mpiilo, mpiihi);
#endif
   }


/*! Set all the elements of a global data structure to zero.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ZERO
 */
   int PPIDD_Zero(int64_t *handle) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_zero(mpihandle);
      else {
         int mproc=NProcs_Work();
         mpierr=twosided_helpga_col(mproc, mpihandle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Zero: array %d has been set to zero.\n",ProcID(),mpihandle);
      if (mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int ihandle = (int) *handle;
      GA_Zero(ihandle);
      return 1 ;
#else
      printf(" ERROR: PPIDD_Zero should not be called in serial case.\n");
      exit(1);
#endif
   }

/*! Get the next shared counter number(helper process should be enabled).
 *
 *  Increment a counter by 1 and returns the counter value (0, 1, ...). */
#ifdef GA_MPI
   static int PPIDD_Nxtval_initialised=0;
   static int64_t PPIDD_Nxtval_handle;
#endif
   void PPIDD_Nxtval(int64_t *numproc, int64_t *val) {
#ifdef MPI2
      if (use_helper_server==0) {
        fprintf(stderr,"%5d: ERROR: Attemp to call NXTVAL routine without helper process!\n", ProcID());
        MPI_Abort(mpigv(Compute_comm),911);
      }
      else {
        int mproc = (int) *numproc;
        *val= (int64_t) NXTVAL(&mproc);
      }
#elif defined(GA_MPI)
      if (*numproc < 0) {
        /* reset - collective */
        if (PPIDD_Nxtval_initialised) PPIDD_Destroy(&PPIDD_Nxtval_handle);
        PPIDD_Nxtval_initialised=0;
        //      }
        //else if (! PPIDD_Nxtval_initialised) {
        /* first call needs to be collective and will return 0*/
        int64_t lentot=1, storetype=1;
        PPIDD_Create(strdup("Nxtval"),&lentot,0,&storetype,&PPIDD_Nxtval_handle);
        PPIDD_Zero(&PPIDD_Nxtval_handle);
        PPIDD_Nxtval_initialised=1;
        *val=0;
      }
      else {
        int64_t inum=1,incr=1;
        PPIDD_Read_inc(&PPIDD_Nxtval_handle,&inum,&incr,val);
      }
#else
      printf(" ERROR: PPIDD_Nxtval should not be called in serial case.\n");
      exit(1);
#endif
   }



/*! Create a new global array by applying all the properties of another existing global.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DUPLICATE
 *  - \b MPI2  does nothing
 */
   void PPIDD_Duplicate(int64_t *handlei, int64_t *handlej, char *name) {
#ifdef GA_MPI
      int ga_a=(int)*handlei;
      int ga_b = GA_Duplicate(ga_a, name);
      *handlej=(int64_t)ga_b;
#else
      printf(" ERROR: PPIDD_Duplicate should not be called in serial and MPI2 cases.\n");
      exit(1);
#endif
}


/*! Returns the name of a global array represented by the handle.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#INQUIRE_NAME
 *  - \c This operation is local.
 */
   void PPIDD_Inquire_name(int64_t *handle, char *name) {
#ifdef MPI2
      char *name2;
      int lxi;
      int i,len_actual;
      int mpihandle = (int) *handle;

      lxi=strlen(name);
      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpiga_inquire_name(mpihandle, &name2);
      else {
         twosided_helpga_inquire_name(mpihandle, &name2);
      }
      len_actual=strlen(name2);
      strncpy(name,name2,len_actual);
      for(i=len_actual;i<lxi;i++) name[i]=' ';
      if(MPIGA_Debug)printf("In PPIDD_Inquire_name: name2=%s,strlen(name2)=%d,lxi=%d\n",name2,len_actual,lxi);
#elif defined(GA_MPI)
      char *name2;
      int lxi;
      int i,len_actual;
      int gahandle = (int) *handle;

      lxi=strlen(name);
/*      strcpy(name2=malloc(80*sizeof(char)),GA_Inquire_name(gahandle));
      strncpy(name,name2,strlen(name2)); */
      name2=GA_Inquire_name(gahandle);
      len_actual=strlen(name2);
      strncpy(name,name2,len_actual);
      for(i=len_actual;i<lxi;i++) name[i]=' ';
      if(MPIGA_Debug)printf("In PPIDD_Inquire_name: name2=%s,strlen(name2)=%d,lxi=%d\n",name2,len_actual,lxi);
#else
      printf(" ERROR: PPIDD_Inquire_name should not be called in serial case.\n");
      exit(1);
#endif
}


/*! Returns the storetype of a global array represented by the handle.
 *
 *     storetype: number of helper processes for storing this global array.
 *  - \c storetype=0 : Normal distributed array stored across the distributed processes
 *  - \c storetype>=1: Low-latency array stored on one or more helper processes (effective only when helper process is enabled).
 *  - \c This operation is local.
 */
   void PPIDD_Inquire_stype(int64_t *handle, int64_t *storetype) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      *storetype=(int64_t)mpiga_inquire_storetype(mpihandle);
#elif defined(GA_MPI)
      *storetype=(int64_t)0;
#else
      printf(" ERROR: PPIDD_Inquire_stype should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Get the amount of memory (in bytes) used in the allocated distributed arrays on the calling processor.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#INQUIRE_NAME
 */
   void PPIDD_Inquire_mem(int64_t *mem_used) {
#ifdef MPI2
      long localmem;
      localmem=mpiga_localmem();
      *mem_used=(int64_t)localmem;
#elif defined(GA_MPI)
      size_t localmem;
      localmem=GA_Inquire_memory();
      *mem_used=(int64_t)localmem;
#else
      *mem_used=(int64_t)0;
#endif
   }



/*! Create a set of mutex variables that can be used for mutual exclusion.
 *
 *  Returns <tt>.true.</tt> if the operation succeeded or <tt>.false.</tt> when failed. It is a collective operation.

 *  - For \b GA, \c storetype doesn't take effect, analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#CREATE_MUTEXES
 *  - For \b MPI2, there are two kinds of mutexes:
 *   - (1) \c storetype=0: mutex data will be stored by a global array across the distributed processes.
 *        Due to the poor performance of MPI-2 one-sided communication, this kind of mutex is very slow.
 *   - (2) \c storetype=1: mutex data will be stored on the helper process. Improve the performance of (1), and all communications are based on two-sided opeartions.
 *
 *   When helper process is disabled, \c storetype doesn't take effect, and mutex data are always stored across the distributed processes.

 */
   int PPIDD_Create_mutexes(int64_t *storetype,int64_t *number) {
#ifdef MPI2
      int stype     = (int) *storetype;
      int mpinumber = (int) *number;
      int mpierr;

      if (use_helper_server==0) {
        mpierr=mpiga_create_mutexes(mpinumber);      /* mutexes data store by a global array across the distributed processes */
      }
      else {
        if (stype==0)
         mpierr=mpiga_create_mutexes(mpinumber);      /* mutexes data store by a global array across the distributed processes */
        else
         mpierr=alloc_general_helpmutexes(mpinumber);  /* mutexes data on helper process */
      }

      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int mpinumber = (int) *number;
      int mpierr;
      mpierr=GA_Create_mutexes(mpinumber);
/* This is one of exceptions in GA (see global/src/capi.c) : Returns [1] if the operation succeeded or [0] when failed */
      if(MPIGA_Debug)printf("In PPIDD_Create_Mutexes: mpierr=%d.\n",mpierr);
      if(mpierr==1) return 1 ;
      else return 0 ;
#else
      return 1 ;
#endif
   }


/*! Lock a mutex object identified by a given mutex number.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#LOCK
 */
   void PPIDD_Lock_mutex(int64_t *inum) {
#ifdef MPI2
      int mpiinum = (int) *inum;
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_lock_mutex(mpiinum);   /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=lock_general_helpmutex(mpiinum); /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpigv(Compute_comm),911);
#elif defined(GA_MPI)
      int mpiinum = (int) *inum;
      GA_Lock(mpiinum);
#endif
   }


/*! Unlock  a mutex object identified by a given mutex number.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#UNLOCK
 */
   void PPIDD_Unlock_mutex(int64_t *inum) {
#ifdef MPI2
      int mpiinum = (int) *inum;
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_unlock_mutex(mpiinum);    /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=unlock_general_helpmutex(mpiinum);  /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpigv(Compute_comm),911);
#elif defined(GA_MPI)
      int mpiinum = (int) *inum;
      GA_Unlock(mpiinum);
#endif
   }


/*! Destroy the set of mutexes created with PPIDD_Create_mutexes.
 *
 * Returns <tt>.true.</tt> if the operation succeeded or <tt>.false.</tt> when failed. This is a collective operation.
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DESTROY_MUTEXES
 */
   int PPIDD_Destroy_mutexes() {
#ifdef MPI2
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_destroy_mutexes();  /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=free_general_helpmutexes();   /* mutexes data on helper process */
      if(mpierr==0) return 1 ;
      else return 0 ;
#elif defined(GA_MPI)
      int mpierr;
      mpierr=GA_Destroy_mutexes();
/* This is one of exceptions in GA (see global/src/capi.c) : Returns [1] if the operation succeeded or [0] when failed */
      if(MPIGA_Debug)printf("In PPIDD_Destroy_Mutexes: mpierr=%d.\n",mpierr);
      if(mpierr==1) return 1 ;
      else return 0 ;
#else
      return 1 ;
#endif
   }

}
