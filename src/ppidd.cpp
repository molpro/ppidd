#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MPI_H
#include <string>
#include <mpi.h>
#endif
#include "ppidd.h"
#include "ppidd_ga_mpi.h"
#include "ppidd_mpi2.h"
#include "ppidd_no_mpi.h"

/*! \file
 * \brief This file contains the PPIDD functions */

static const int ppidd_impl_default= PPIDD_IMPL_DEFAULT;
static const int ppidd_impl_no_mpi = PPIDD_IMPL_NO_MPI;
static const int ppidd_impl_ga_mpi = PPIDD_IMPL_GA_MPI;
static const int ppidd_impl_mpi2   = PPIDD_IMPL_MPI2;

static int ppidd_impl=ppidd_impl_default;

extern "C" {

/*! \brief Initialize the PPIDD parallel environment
    \details
    - For \b GA, includes initialization of MPI and GA.
    - For \b MPI2, calls MPI_Init. */
   void PPIDD_Initialize(int *argc, char ***argv, int impl) {
    switch (impl) {
     case (ppidd_impl_ga_mpi):
     case (ppidd_impl_mpi2):
     case (ppidd_impl_no_mpi):
      break;
     default:
      fprintf(stderr,"ERROR: impl '%d' is unknown\n",impl);
      exit(1);
    }

    switch (impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
#endif
     case (ppidd_impl_mpi2):
#endif
     case (ppidd_impl_no_mpi):
      ppidd_impl=impl;
      break;
     default:
      fprintf(stderr,"ERROR: impl '%d' is unavailable\n",impl);
      exit(1);
    }

#ifdef HAVE_MPI_H
    int flag=0;
    int ret=MPI_Initialized(&flag);
    if (ret != MPI_SUCCESS) {fprintf(stderr,"MPI_Initialized failed (%d)",ret); exit(1);}
    if (flag) {std::string msg="MPI already initialized"; PPIDD_Error(&msg[0],flag);}
    ret=MPI_Init(argc, argv);
    if (ret != MPI_SUCCESS) {fprintf(stderr,"MPI_Init failed (%d)",ret); exit(1);}
#ifdef HAVE_GA_H
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) ppidd_impl=ppidd_impl_mpi2; /* for single process switch to mpi2 version (because otherwise GA built with mpi-pr would fail) */
#endif
#endif

    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Initialize(argc,argv,impl);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Initialize(argc,argv,impl);
#endif
     default:
      return no_mpi::PPIDD_Initialize(argc,argv,impl);
    }
   }


/*! Initialize the PPIDD data structure
 *
 *  - For \b GA, does nothing
 *  - For \b MPI2, Initialize global data structure and set helper server.
 */
   void PPIDD_Initialize_data() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Initialize_data();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Initialize_data();
#endif
     default:
      return no_mpi::PPIDD_Initialize_data();
    }
   }


/*! Return communicator for worker process group, which excludes any helper process.
 * For simplicity, the returned data type is set as integer. So for C calling, it may have to be converted
 * into MPI_Comm by MPI_Comm_f2c afterward if needed.
 *
 *  - For \b MPI2, Return communicator for worker process group.
 *  - For \b GA and serial cases, should not be called.
 */
   int64_t PPIDD_Worker_comm() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Worker_comm();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Worker_comm();
#endif
     default:
      return no_mpi::PPIDD_Worker_comm();
    }
   }


/*! Terminate the PPIDD parallel environment.
 *
 *  - For \b GA, tidy up global arrays and MPI, analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#TERMINATE
 *  - For \b MPI2, tidy up some associated resources and call MPI_Finalize.
 */
   void PPIDD_Finalize() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Finalize();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Finalize();
#endif
     default:
      return no_mpi::PPIDD_Finalize();
    }
   }


/*  Detect whether MA is used for allocation of GA memory.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#USES_MA
 *  - \b MPI2 always returns <tt>.false.</tt>
 */
   int PPIDD_Uses_ma() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Uses_ma();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Uses_ma();
#endif
     default:
      return no_mpi::PPIDD_Uses_ma();
    }
   }


/*  Initialize the memory allocator.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/ma/MAapi.html
 *  - \b MPI2 always returns <tt>.true.</tt>
 */
   int PPIDD_MA_init(int dtype, int64_t *stack, int64_t *heap) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_MA_init(dtype,stack,heap);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_MA_init(dtype,stack,heap);
#endif
     default:
      return no_mpi::PPIDD_MA_init(dtype,stack,heap);
    }
   }


/*! Return an elapsed time on the calling process.
 *
 * Returns a floating-point number of seconds, representing elapsed wall-clock time since an arbitrary time in the past.
 * The times returned are local to the node that called them. There is no requirement that different nodes return the same time.
 * This is a local operation.
 *
 *  - \b GA calls GA_Wtime, which is available only in release 4.1 or greater, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#WTIME
 *  - \b MPI2 calls MPI_Wtime
 */
   double PPIDD_Wtime() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Wtime();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Wtime();
#endif
     default:
      return no_mpi::PPIDD_Wtime();
    }
  }


/*! Print an error message and abort the program execution.
 *
 *  - \b GA calls GA_Error, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ERROR
 *  - For \b MPI2, prints error, and then calls MPI_Abort.
 */
   void PPIDD_Error(char *message, int code) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Error(message,code);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Error(message,code);
#endif
     default:
      return no_mpi::PPIDD_Error(message,code);
    }
   }


/*! Set the flag for data helper server and set the number of how many processes own one helper server
 *
 *  Set the helper_server flag: 1 (use); 0 (don't use). */
   void PPIDD_Helper_server(int flag, int numprocs_per_server) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Helper_server(flag,numprocs_per_server);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Helper_server(flag,numprocs_per_server);
#endif
     default:
      return no_mpi::PPIDD_Helper_server(flag,numprocs_per_server);
    }
   }


/*!  Determine the total number of processes available (including helper process if there is one).
 *
 *  - \b GA calls GA_Nnodes,
 *  http://hpc.pnl.gov/globalarrays/api/c_op_api.html#NNODES
 *  - \b MPI2 calls MPI_Comm_size for communicator MPI_COMM_WORLD.
 */
   int PPIDD_Size_all() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Size_all();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Size_all();
#endif
     default:
      return no_mpi::PPIDD_Size_all();
    }
   }


/*  Determine the number of compute processes.
 *
 *  - \b GA calls GA_Nnodes,
 *  http://hpc.pnl.gov/globalarrays/api/c_op_api.html#NNODES
 *  - \b MPI2 calls MPI_Comm_size for computational communicator.
 */
   int PPIDD_Size() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Size();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Size();
#endif
     default:
      return no_mpi::PPIDD_Size();
    }
   }


/*! Determine the rank of the calling process.
 *
 *  - \b GA calls GA_Nodeid,
 *  http://hpc.pnl.gov/globalarrays/api/c_op_api.html#NODEID
 *  - \b MPI2 calls MPI_Comm_rank in computational communicator.
 */
   int PPIDD_Rank() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Rank();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Rank();
#endif
     default:
      return no_mpi::PPIDD_Rank();
    }
   }


/*! Initialize tracing of completion status of data movement operations.
 *
 *  - \b GA calls GA_Init_fence, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#INIT_FENCE
 *  - \b MPI2 does nothing
 */
   void PPIDD_Init_fence() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Init_fence();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Init_fence();
#endif
     default:
      return no_mpi::PPIDD_Init_fence();
    }
   }


/*! Block the calling process until all data transfers complete.
 *
 *  - \b GA calls GA_Fence, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#FENCE
 *  - \b MPI2 does nothing
 */
   void PPIDD_Fence() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Fence();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Fence();
#endif
     default:
      return no_mpi::PPIDD_Fence();
    }
   }


/*! Blocking/nonblocking send.
 *
 *  - \b GA calls SND_.
 *  - \b MPI2 calls MPI_Send ( sync is 1) or MPI_Isend ( sync is 0).
 */
   void PPIDD_Send(void *buf,int64_t *count,int dtype,int64_t *dest,int sync) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Send(buf,count,dtype,dest,sync);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Send(buf,count,dtype,dest,sync);
#endif
     default:
      return no_mpi::PPIDD_Send(buf,count,dtype,dest,sync);
    }
   }


/*! Blocking/nonblocking receive.
 *
 *  - \b GA calls RCV_.
 *  - \b MPI2 calls MPI_Recv ( sync is 1) or MPI_Irecv( sync is 0).
 */
   void PPIDD_Recv(void *buf,int64_t *count,int dtype,int64_t *source,int64_t *lenreal,int64_t *sourcereal,int sync) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Recv(buf,count,dtype,source,lenreal,sourcereal,sync);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Recv(buf,count,dtype,source,lenreal,sourcereal,sync);
#endif
     default:
      return no_mpi::PPIDD_Recv(buf,count,dtype,source,lenreal,sourcereal,sync);
    }
   }


/*! Wait for completion of all asynchronous send/receive.
 *
 *  - \b MPI2/GA_MPI calls MPI_Wait for all asynchronous requests.
 */
   void PPIDD_Wait() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Wait();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Wait();
#endif
     default:
      return no_mpi::PPIDD_Wait();
    }
   }


/*! Detect whether a message of the specified type is available.
 *
 *  Return <tt>.true.</tt> if the message is available, otherwise <tt>.false.</tt>.
 *
 *  - \b MPI2/GA_MPI calls MPI_Iprobe
 */
   int PPIDD_Iprobe(int tag,int64_t *source) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Iprobe(tag,source);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Iprobe(tag,source);
#endif
     default:
      return no_mpi::PPIDD_Iprobe(tag,source);
    }
   }


/*! Broadcast a message from the root process to all other processes.
 *
 *  Collective operation.
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#BRDCST
 *  - \b MPI2 calls MPI_Bcast
 *
 *  - \c dtype=0 : Fortran integer and logical types
 *  - \c dtype=1 : Fortran double precision type */
   void PPIDD_BCast(void *buffer,int64_t *count,int dtype,int root) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_BCast(buffer,count,dtype,root);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_BCast(buffer,count,dtype,root);
#endif
     default:
      return no_mpi::PPIDD_BCast(buffer,count,dtype,root);
    }
   }


/*! Synchronize processes and ensure all have reached this routine.
 *
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#SYNC
 *  - \b MPI2 calls MPI_Barrier
 */
   void PPIDD_Barrier() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Barrier();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Barrier();
#endif
     default:
      return no_mpi::PPIDD_Barrier();
    }
   }


/*! Combine values from all processes and distribute the result back to all processes.
 *
 *  - \b GA analogous to GA_dgop and GA_igop, http://hpc.pnl.gov/globalarrays/api/c_op_api.html#GOP
 *  - For \b MPI2, it is realized by MPI_Allreduce
 *
 *  - \c type=0 : Fortran Integer
 *  - \c type=1 : Fortran Double Precision */
   void PPIDD_Gsum(int dtype,void *buffer,int len, char *op) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Gsum(dtype,buffer,len,op);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Gsum(dtype,buffer,len,op);
#endif
     default:
      return no_mpi::PPIDD_Gsum(dtype,buffer,len,op);
    }
   }


/*! Create an array by following the user-specified distribution and return integer handle representing the array.
 *
 * Irregular distributed array data are stored across the distributed processes.
 *  - \c dtype=0 : Fortran integer and logical types
 *  - \c dtype=1 : Fortran double precision type

 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#CREATE_IRREG
 */
   int PPIDD_Create_irreg(char *name, int64_t *lenin, int64_t nchunk, int dtype, int storetype, int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Create_irreg(name,lenin,nchunk,dtype,storetype,handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Create_irreg(name,lenin,nchunk,dtype,storetype,handle);
#endif
     default:
      return no_mpi::PPIDD_Create_irreg(name,lenin,nchunk,dtype,storetype,handle);
    }
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
   int PPIDD_Create(char *name,int64_t lentot, int dtype, int storetype, int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Create(name,lentot,dtype,storetype,handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Create(name,lentot,dtype,storetype,handle);
#endif
     default:
      return no_mpi::PPIDD_Create(name,lentot,dtype,storetype,handle);
    }
   }


/*! Deallocate the array represented by handle and free any associated resources.
 *
 *  - \b GA analogous http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DESTROY
 */
   int PPIDD_Destroy(int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Destroy(handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Destroy(handle);
#endif
     default:
      return no_mpi::PPIDD_Destroy(handle);
    }
   }


/*! Return the range of a distributed array held by a specified process.
 *
 *  Return <tt>.true.</tt> if successful, otherwise <tt>.false.</tt>
 *  If no array elements are owned by the process, the range is returned as [0,-1].
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DISTRIBUTION
 */
   int PPIDD_Distrib(int64_t *handle,int rank,int64_t *ilo,int64_t *ihi) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Distrib(handle,rank,ilo,ihi);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Distrib(handle,rank,ilo,ihi);
#endif
     default:
      return no_mpi::PPIDD_Distrib(handle,rank,ilo,ihi);
    }
   }


/*! Return a list of the processes that hold the data.
 *
 *  Parts of the specified patch might be actually 'owned' by several processes.
 *  np is the number of processes hold tha data (return 0  if ilo/ihi are out of bounds "0").
 *  - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#LOCATE_REGION
 */
   int PPIDD_Location(int64_t *handle,     /*!< array handle */
                       int64_t *ilo,        /*!< lower element subscript, 1 (not 0) for first element */
                       int64_t *ihi,        /*!< higher element subscript */
                       int64_t *map,        /*!< return start/end index for <tt>proclist</tt> */
                       int64_t *proclist,   /*!< proc id list */
                       int64_t *np          /*!< proc number */
   ) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Location(handle,ilo,ihi,map,proclist,np);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Location(handle,ilo,ihi,map,proclist,np);
#endif
     default:
      return no_mpi::PPIDD_Location(handle,ilo,ihi,map,proclist,np);
    }
   }


/*! \brief Copies data from array section to the local array buffer according to starting and ending index.
    \details analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#GET */
   int PPIDD_Get(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Get(handle,ilo,ihi,buff);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Get(handle,ilo,ihi,buff);
#endif
     default:
      return no_mpi::PPIDD_Get(handle,ilo,ihi,buff);
    }
   }

/*! \brief Put local buffer data into a section of a global array according to starting and ending index.
    \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#PUT */
   int PPIDD_Put(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Put(handle,ilo,ihi,buff);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Put(handle,ilo,ihi,buff);
#endif
     default:
      return no_mpi::PPIDD_Put(handle,ilo,ihi,buff);
    }
   }


/*! \brief Accumulate data into a section of a global array.
    \details Atomic operation. global array section (ilo, ihi) += *fac * buffer.
    Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ACC */
   int PPIDD_Acc(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Acc(handle,ilo,ihi,buff,fac);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Acc(handle,ilo,ihi,buff,fac);
#endif
     default:
      return no_mpi::PPIDD_Acc(handle,ilo,ihi,buff,fac);
    }
   }


/*! \brief Atomically read and increment an element in an integer array.
    \details Reads data from the (inum) element of a global array of integers, returns that value, and increments the (inum)
    element by a given increment. This is a fetch-and-add operation.
    Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#READ_INC */
   void PPIDD_Read_inc(int64_t *handle,int64_t *inum,int64_t *incr,int64_t *returnval) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Read_inc(handle,inum,incr,returnval);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Read_inc(handle,inum,incr,returnval);
#endif
     default:
      return no_mpi::PPIDD_Read_inc(handle,inum,incr,returnval);
    }
   }


/*! \brief Set all the elements in an array patch to zero.
    \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ZERO_PATCH */
   void PPIDD_Zero_patch(int64_t *handle,int64_t *ilo,int64_t *ihi) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Zero_patch(handle,ilo,ihi);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Zero_patch(handle,ilo,ihi);
#endif
     default:
      return no_mpi::PPIDD_Zero_patch(handle,ilo,ihi);
    }
   }


/*! \brief Set all the elements of a global data structure to zero.
 *  \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#ZERO */
   int PPIDD_Zero(int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Zero(handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Zero(handle);
#endif
     default:
      return no_mpi::PPIDD_Zero(handle);
    }
   }

/*! \brief Get the next shared counter number(helper process should be enabled).
    \details Increment a counter by 1 and returns the counter value (0, 1, ...). */
   void PPIDD_Nxtval(int numproc, int64_t *val) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Nxtval(numproc,val);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Nxtval(numproc,val);
#endif
     default:
      return no_mpi::PPIDD_Nxtval(numproc,val);
    }
   }


/*! \brief Create a new global array by applying all the properties of another existing global.
    \details
    - \b GA analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DUPLICATE
    - \b MPI2 does nothing */
   void PPIDD_Duplicate(int64_t *handlei, int64_t *handlej, char *name) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Duplicate(handlei,handlej,name);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Duplicate(handlei,handlej,name);
#endif
     default:
      return no_mpi::PPIDD_Duplicate(handlei,handlej,name);
    }
   }


/*! \brief Returns the name of a global array represented by the handle.
    \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#INQUIRE_NAME
    This operation is local. */
   void PPIDD_Inquire_name(int64_t *handle, char *name) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Inquire_name(handle,name);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Inquire_name(handle,name);
#endif
     default:
      return no_mpi::PPIDD_Inquire_name(handle,name);
    }
   }


/* \brief Returns the storetype of a global array represented by the handle.
   \details storetype: number of helper processes for storing this global array.
   - \c storetype=0 : Normal distributed array stored across the distributed processes
   - \c storetype>=1: Low-latency array stored on one or more helper processes (effective only when helper process is enabled).
   - \c This operation is local. */
   void PPIDD_Inquire_stype(int64_t *handle, int64_t *storetype) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Inquire_stype(handle,storetype);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Inquire_stype(handle,storetype);
#endif
     default:
      return no_mpi::PPIDD_Inquire_stype(handle,storetype);
    }
   }


/*! \brief Get the amount of memory (in bytes) used in the allocated distributed arrays on the calling processor.
    \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#INQUIRE_NAME */
   void PPIDD_Inquire_mem(int64_t *mem_used) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Inquire_mem(mem_used);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Inquire_mem(mem_used);
#endif
     default:
      return no_mpi::PPIDD_Inquire_mem(mem_used);
    }
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
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Create_mutexes(storetype,number);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Create_mutexes(storetype,number);
#endif
     default:
      return no_mpi::PPIDD_Create_mutexes(storetype,number);
    }
   }


/*! \brief Lock a mutex object identified by a given mutex number.
    \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#LOCK */
   void PPIDD_Lock_mutex(int64_t *inum) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Lock_mutex(inum);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Lock_mutex(inum);
#endif
     default:
      return no_mpi::PPIDD_Lock_mutex(inum);
    }
   }


/*! \brief Unlock a mutex object identified by a given mutex number.
    \details Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#UNLOCK */
   void PPIDD_Unlock_mutex(int64_t *inum) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Unlock_mutex(inum);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Unlock_mutex(inum);
#endif
     default:
      return no_mpi::PPIDD_Unlock_mutex(inum);
    }
   }


/*! \brief Destroy the set of mutexes created with PPIDD_Create_mutexes.
    \return 1 if the operation succeeded or 0 when failed.
    \details This is a collective operation. Analogous to http://hpc.pnl.gov/globalarrays/api/c_op_api.html#DESTROY_MUTEXES */
   int PPIDD_Destroy_mutexes() {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Destroy_mutexes();
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Destroy_mutexes();
#endif
     default:
      return no_mpi::PPIDD_Destroy_mutexes();
    }
   }


/* ************************************************************************ *\
   Creates an EAF file using name and path specified in name as a template.
   Return the EAF file descriptor in handle.
   It is a non-collective operation.
\* ************************************************************************ */
   int PPIDD_Eaf_open(char *name,int64_t *type, int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_open(name,type,handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_open(name,type,handle);
#endif
     default:
      return no_mpi::PPIDD_Eaf_open(name,type,handle);
    }
   }


/* ******************************************************************************************** *\
   Synchronously write to the file specified by the file handle.
   Writes number of bytes to the file identified by handle at location offset.
\* ******************************************************************************************** */
   int PPIDD_Eaf_write(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_write(handle,byte_offset,buff,byte_length);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_write(handle,byte_offset,buff,byte_length);
#endif
     default:
      return no_mpi::PPIDD_Eaf_write(handle,byte_offset,buff,byte_length);
    }
   }


/* ******************************************************************************************** *\
   Asynchronously write to the file specified by the file handle, and return a handle to the asynchronous operation.
   Writes number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when eaf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   int PPIDD_Eaf_awrite(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_awrite(handle,byte_offset,buff,byte_length,request_id);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_awrite(handle,byte_offset,buff,byte_length,request_id);
#endif
     default:
      return no_mpi::PPIDD_Eaf_awrite(handle,byte_offset,buff,byte_length,request_id);
    }
   }


/* ******************************************************************************************** *\
   Synchronously read from the file specified by the file handle.
   Reads number of bytes to the file identified by handle at location offset.
\* ******************************************************************************************** */
   int PPIDD_Eaf_read(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_read(handle,byte_offset,buff,byte_length);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_read(handle,byte_offset,buff,byte_length);
#endif
     default:
      return no_mpi::PPIDD_Eaf_read(handle,byte_offset,buff,byte_length);
    }
   }


/* ******************************************************************************************** *\
   Asynchronously read from the file specified by the file handle, and return a handle to the asynchronous operation.
   Reads number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when eaf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   int PPIDD_Eaf_aread(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_aread(handle,byte_offset,buff,byte_length,request_id);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_aread(handle,byte_offset,buff,byte_length,request_id);
#endif
     default:
      return no_mpi::PPIDD_Eaf_aread(handle,byte_offset,buff,byte_length,request_id);
    }
   }


/* ************************************************************************************ *\
   Wait for the completion of the asynchronous request associated with request_id.
   Request_id is invalidated after the calling.
   integer request_id   --[in]  Handle of asynchronous request.
   integer ierr         --[out] Error code. 0 if it is able to wait for completion,
                          else returns error code.
\* ************************************************************************************ */
   int PPIDD_Eaf_wait(int64_t *handle,int64_t *request_id) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_wait(handle,request_id);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_wait(handle,request_id);
#endif
     default:
      return no_mpi::PPIDD_Eaf_wait(handle,request_id);
    }
   }


/* ********************************************************************************** *\
   Blocks the calling process until all of the num I/O operations associated with ids
   specified in list complete. Finally invalidates (modifies) ids on the list.
\* ********************************************************************************** */
   int PPIDD_Eaf_waitall(int64_t *list, int64_t *num) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_waitall(list,num);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_waitall(list,num);
#endif
     default:
      return no_mpi::PPIDD_Eaf_waitall(list,num);
    }
   }


/* ************************************************************************************ *\
   Determine if an asynchronous request has completed or is pending.
   integer request_id   --[in]  Handle of asynchronous request.
   integer status       --[out] Pending or completed status argument.
                          status returns 0 if the asynchronous operation is complete, or 1 otherwise.
                          If the asynchronous request is complete, id is invalidated.
   integer ierr         --[out] Error code. 0 if probe succeeded, else returns error code.
\* ************************************************************************************ */
   int PPIDD_Eaf_probe(int64_t *request_id,int64_t *status) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_probe(request_id,status);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_probe(request_id,status);
#endif
     default:
      return no_mpi::PPIDD_Eaf_probe(request_id,status);
    }
   }


/* ************************************************************************************ *\
   Close an EAF file.
   integer handle  --[in]  File Handle.
   integer ierr    --[out] Error code. 0 if the file was closed, else returns error code.
\* ************************************************************************************ */
   int PPIDD_Eaf_close(int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_close(handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_close(handle);
#endif
     default:
      return no_mpi::PPIDD_Eaf_close(handle);
    }
   }


/* ********************************************************************************************* *\
   Delete an EAF file
   character*(*) name   -- [in]  File name.
   integer ierr         -- [out] Error code. 0 if the file was deleted, else returns error code.
\* ********************************************************************************************* */
   int PPIDD_Eaf_delete(char *name) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_delete(name);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_delete(name);
#endif
     default:
      return no_mpi::PPIDD_Eaf_delete(name);
    }
   }


/* ************************************************************************************ *\
   Determine the length (in bytes) of an EAF file.
   integer handle    --[in]  File Handle.
   double fsize      --[out] File length in bytes.
\* ************************************************************************************ */
   int PPIDD_Eaf_length(int64_t *handle,double *fsize) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_length(handle,fsize);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_length(handle,fsize);
#endif
     default:
      return no_mpi::PPIDD_Eaf_length(handle,fsize);
    }
   }


/* *************************************************************************************** *\
   Truncate an EAF file at specified offset (in bytes).
   integer handle --[in]  File Handle.
   double offset  --[in]  Offset in bytes.
   integer ierr   --[out] Error code. 0 if the file was truncated, else returns error code.
\* *************************************************************************************** */
   int PPIDD_Eaf_truncate(int64_t *handle,double *offset) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_truncate(handle,offset);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_truncate(handle,offset);
#endif
     default:
      return no_mpi::PPIDD_Eaf_truncate(handle,offset);
    }
   }


/* ********************************************************************************* *\
   Returns a string interpretation of the error code, or an empty string
   (Fortran all blanks, C null first character) if the error code is not recognized.
        code             -- [in]  error code returned by a previous call to EAF
        message          -- [out] character string where the corresponding message
\* ********************************************************************************* */
   void PPIDD_Eaf_errmsg(int *code,char *message) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Eaf_errmsg(code,message);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Eaf_errmsg(code,message);
#endif
     default:
      return no_mpi::PPIDD_Eaf_errmsg(code,message);
    }
   }


/* ************************************************************************ *\
   Creates shared file using name and path specified in name as a template.
   req_size specifies size of a typical request (-1. means "don't know").
   It is a collective operation.
\* ************************************************************************ */
   int PPIDD_Sf_create(char *name ,double *size_hard_limit, double *size_soft_limit, double *req_size, int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_create(name,size_hard_limit,size_soft_limit,req_size,handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_create(name,size_hard_limit,size_soft_limit,req_size,handle);
#endif
     default:
      return no_mpi::PPIDD_Sf_create(name,size_hard_limit,size_soft_limit,req_size,handle);
    }
   }


/* ******************************************************************************************** *\
   Asynchronous write operation.
   Writes number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when sf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   int PPIDD_Sf_write(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_write(handle,byte_offset,byte_length,buff,request_id);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_write(handle,byte_offset,byte_length,buff,request_id);
#endif
     default:
      return no_mpi::PPIDD_Sf_write(handle,byte_offset,byte_length,buff,request_id);
    }
   }

/* ******************************************************************************************** *\
   Asynchronous read operation.
   Reads number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when sf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   int PPIDD_Sf_read(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_read(handle,byte_offset,byte_length,buff,request_id);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_read(handle,byte_offset,byte_length,buff,request_id);
#endif
     default:
      return no_mpi::PPIDD_Sf_read(handle,byte_offset,byte_length,buff,request_id);
    }
   }


/* ************************************************************************************ *\
   Blocks the calling process until I/O operation associated with request_id completes.
   Invalidates request_id.
\* ************************************************************************************ */
   int PPIDD_Sf_wait(int64_t *request_id) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_wait(request_id);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_wait(request_id);
#endif
     default:
      return no_mpi::PPIDD_Sf_wait(request_id);
    }
   }


/* ********************************************************************************** *\
   Blocks the calling process until all of the num I/O operations associated with ids
   specified in list complete. Invalidates (modifies) ids on the list.
\* ********************************************************************************** */
   int PPIDD_Sf_waitall(int64_t *list, int64_t *num) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_waitall(list,num);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_waitall(list,num);
#endif
     default:
      return no_mpi::PPIDD_Sf_waitall(list,num);
    }
   }


/* ************************************************ *\
   Destroys the shared file associated with handle.
   It is a collective operation.
\* ************************************************ */
   int PPIDD_Sf_destroy(int64_t *handle) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_destroy(handle);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_destroy(handle);
#endif
     default:
      return no_mpi::PPIDD_Sf_destroy(handle);
    }
   }


/* ********************************************************************************* *\
   Returns a string interpretation of the error code, or an empty string
   (Fortran all blanks, C null first character) if the error code is not recognized.
        code             -- error code returned by a previous call to SF [in]
        message          -- character string where the corresponding message
                            is written [out]
\* ********************************************************************************* */
   void PPIDD_Sf_errmsg(int *code,char *message) {
    switch (ppidd_impl) {
#ifdef HAVE_MPI_H
#ifdef HAVE_GA_H
     case (ppidd_impl_ga_mpi):
      return ga_mpi::PPIDD_Sf_errmsg(code,message);
#endif
     case (ppidd_impl_mpi2):
      return mpi2::PPIDD_Sf_errmsg(code,message);
#endif
     default:
      return no_mpi::PPIDD_Sf_errmsg(code,message);
    }
   }

}
