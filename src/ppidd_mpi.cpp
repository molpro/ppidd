#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif
#ifdef HAVE_MPI_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <climits>
#include <vector>
#include "mpiga_base.h"
#include "mpi_nxtval.h"
#include "mpi_utils.h"
#include "ppidd_mpi.h"
#define   MPI_EAF_RW -1
#define   MPI_EAF_W  -2
#define   MPI_EAF_R  -3

namespace ppidd {
namespace mpi {

   static int MPIGA_Debug=0;
   static int MPI_Debug=0;

   [[ noreturn ]] static void do_not_call(const char* function) {
    fprintf(stderr,"%s should not be called in mpi case\n",function);
    exit(1);
   }

   void PPIDD_Initialize(int *argc, char ***argv, int impl, int fortint_size) {
    int mpierr=mpiga_initialize(argc,argv);
    mpi_test_status("PPIDD_Initialize:",mpierr);
    if (fortint_size==4) dtype_mpi[PPIDD_FORTINT] = MPI_INT32_T;
    else if (fortint_size==8) dtype_mpi[PPIDD_FORTINT] = MPI_INT64_T;
    else {
     MPIGA_Error(" PPIDD_Initialize: unable to map PPIDD_FORTINT ", fortint_size);
    }
    dtype_mpi[PPIDD_DOUBLE] = MPI_DOUBLE;
    dtype_mpi[PPIDD_INT] = MPI_INT;
   }


   void PPIDD_Initialize_data() {
      mpiga_initialize_data();
   }


   int PPIDD_Worker_comm() {
      MPI_Comm mycomm=mpiga_compute_comm();

/* test whether worker communicator contains all the processes, if so then return MPI_COMM_WORLD */
      int np_all, np_worker=mpigv::nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &np_all);
      if(np_all==np_worker) mycomm=MPI_COMM_WORLD;

      return MPI_Comm_c2f(mycomm);
   }


   void PPIDD_Finalize() {
      mpiga_terminate();
   }


   int PPIDD_Uses_ma() {
    return 0;
   }


   int PPIDD_MA_init(int dtype, int64_t stack, int64_t heap) {
      return 1;
   }


   double PPIDD_Wtime() {
      return MPI_Wtime();
   }


   void PPIDD_Error(const char *message, int code) {
      MPIGA_Error(message,code);
   }


   void PPIDD_Helper_server(int flag, int numprocs_per_server) {
      if (flag) {                                    /* mutilple helper servers, node helper server, and single helper server */
         use_helper_server=1;
         if ( numprocs_per_server > 1 ) {            /* reasonable mutilple helper servers, and single helper server */
            NPROCS_PER_HELPER=numprocs_per_server;
         }
         else if ( numprocs_per_server == 0 ) {      /* node helper server: one helper server on every node */
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
            fprintf(stderr,"%5d: WARNING: nprocs_per_server=%d is unreasonable. Will use only single helper server for all processes!\n", ProcID(),numprocs_per_server);
            NPROCS_PER_HELPER=99999999;
         }
      }
      else {                                         /* no helper server */
         use_helper_server=0;
         NPROCS_PER_HELPER=-1;
      }
   }


   int PPIDD_Size_all() {
      int mpinp;
      MPI_Comm_size(MPI_COMM_WORLD,&mpinp);
      return mpinp;
   }


   int PPIDD_Size() {
      int mpinp;
      MPI_Comm_size(mpiga_compute_comm(),&mpinp);
      return mpinp;
   }


   int PPIDD_Rank() {
      int mpime;
      MPI_Comm_rank(mpiga_compute_comm(),&mpime);
      return mpime;
   }


   void PPIDD_Init_fence() {
   }


   void PPIDD_Fence() {
   }


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


   void PPIDD_Send(void *buf,int count,int dtype,int dest,int sync) {
      MPI_Comm mpicomm=mpiga_compute_comm();
      int mpitag=dtype;
      int mpierr;

      int mpilenbuf;
      MPI_Type_size( dtype_mpi[dtype], &mpilenbuf);
      mpilenbuf *= count;

      if (MPIGA_Debug) {
         printf("PPIDD_SEND: node %d sending to %d, len(bytes)=%d, mes tag=%d, sync=%d\n",
                 ProcID(), dest, mpilenbuf, mpitag, sync);
         fflush(stdout);
      }

      if (sync) {
         mpierr=MPI_Send(buf,mpilenbuf,MPI_CHAR,dest,mpitag,mpicomm);
         mpi_test_status("PPIDD_SEND: SEND:",mpierr);
      }
      else {
         if (n_in_msg_mpiq >= MAX_MPIQ_LEN) {
            MPIGA_Error("PPIDD_SEND: nonblocking SEND: overflowing async Queue limit",n_in_msg_mpiq);
         }
         mpierr = MPI_Isend(buf, mpilenbuf, MPI_CHAR,dest, mpitag,mpicomm,
                     &msg_mpiq[n_in_msg_mpiq].request);
         mpi_test_status("PPIDD_SEND: nonblocking SEND:",mpierr);

         msg_mpiq[n_in_msg_mpiq].node   =(long) dest;
         msg_mpiq[n_in_msg_mpiq].type   =(long) dtype;
         msg_mpiq[n_in_msg_mpiq].lenbuf =(long) mpilenbuf;
         msg_mpiq[n_in_msg_mpiq].snd = (long)1;
      }
   }


   void PPIDD_Recv(void *buf,int count,int dtype,int source,int *lenreal,int *sourcereal,int sync) {
      MPI_Comm mpicomm=mpiga_compute_comm();
      int mpitag=dtype;
      int mpinode,mpierr;
      MPI_Status status;
      MPI_Request request;

      int mpilenbuf;
      MPI_Type_size( dtype_mpi[dtype], &mpilenbuf);
      mpilenbuf *= count;

      if (source == -1)
         mpinode = MPI_ANY_SOURCE;
      else
         mpinode = source;

      if (MPIGA_Debug) {
         printf("PPIDD_Recv: node %d receving from %d, len(bytes)=%d, mes tag=%d, sync=%d\n",
                 ProcID(), source, mpilenbuf, mpitag, sync);
         fflush(stdout);
      }

      if(sync==0){
         if (n_in_msg_mpiq >= MAX_MPIQ_LEN) {
            MPIGA_Error("PPIDD_Recv: nonblocking RECV: overflowing async Queue limit", n_in_msg_mpiq);
         }
         mpierr = MPI_Irecv(buf,mpilenbuf,MPI_CHAR,mpinode,mpitag,mpicomm,&request);
         mpi_test_status("PPIDD_Recv: nonblocking RECV:",mpierr);

         *sourcereal = mpinode;          /* Get source node  */
         *lenreal = -1;
         msg_mpiq[n_in_msg_mpiq].request = request;
         msg_mpiq[n_in_msg_mpiq].node   = (long)source;
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
         *sourcereal = status.MPI_SOURCE;
         *lenreal    = mpilenbuf;
      }
   }


   void PPIDD_Wait() {
      MPI_Status status;

      for (int i=0; i<n_in_msg_mpiq; i++){
         if (MPIGA_Debug) {
            printf("PPIDD_Wait: node %d waiting for msg to/from %ld, #%d\n", ProcID(), msg_mpiq[i].node, i);
            fflush(stdout);
         }
         int mpierr = MPI_Wait(&msg_mpiq[i].request, &status);
         mpi_test_status("PPIDD_Wait:",mpierr);
      }
      n_in_msg_mpiq = 0;
   }


   int PPIDD_Iprobe(int tag, int source) {
      MPI_Comm mpicomm=mpiga_compute_comm();
      int flag;
      MPI_Status status;

      int mpisource = (source < 0) ? MPI_ANY_SOURCE : source;
      int mpierr = MPI_Iprobe(mpisource, tag, mpicomm, &flag, &status);
      mpi_test_status("PPIDD_Iprobe:",mpierr);
      return flag ? 1 : 0;
   }


   void PPIDD_BCast(void *buffer,int count,int dtype,int root) {
      char *cbuf=(char *)buffer;
      int dtype_size;
      MPI_Type_size( dtype_mpi[dtype], &dtype_size);
      for (int64_t remaining=count, addr=0; remaining > 0; remaining-=(int64_t)BCAST_BATCH_SIZE, addr+=(int64_t)BCAST_BATCH_SIZE) {
       int mpierr=MPI_Bcast(&cbuf[addr*dtype_size],(int)std::min((int64_t)BCAST_BATCH_SIZE,remaining),dtype_mpi[dtype],root,mpiga_compute_comm());
       mpi_test_status("PPIDD_BCast:",mpierr);
      }
   }


   void PPIDD_Barrier() {
      int mpierr=MPI_Barrier(mpiga_compute_comm());
      mpi_test_status("PPIDD_Barrier:",mpierr);
   }


   void PPIDD_Gsum(int dtype,void *buffer,int len, const char *op) {
      MPI_GSum(dtype_mpi[dtype], buffer, len, op);
   }


   int PPIDD_Create_irreg(const char *name, const int64_t *lenin, int64_t nchunk, int dtype, int storetype, int *handle) {
      int mpierr;
      int mpinchunk=(int)nchunk;

      std::vector<int> mpilenin(mpinchunk);
      for (int i=0;i<mpinchunk;i++) mpilenin[i]=(int)lenin[i];
      MPI_Datatype mpidtype=dtype_mpi[dtype];
      if (use_helper_server==0) {
        mpierr=mpiga_create_irreg(mpilenin.data(), mpinchunk, mpidtype, handle);
      }
      else {
        if (storetype==0)
          mpierr=mpiga_create_irreg(mpilenin.data(), mpinchunk, mpidtype, handle);
        else {
          int mproc=0;
          mpierr=twosided_helpga_create_irreg(mproc, mpilenin.data(), mpinchunk, handle, dtype);
        }
      }
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Create(const char *name,int64_t lentot, int dtype, int storetype, int *handle) {
      int mpierr;
      int mpilentot;

      if (lentot > INT_MAX) {
       printf(" ERROR: PPIDD_Create: lentot too large for MPI\n");
       exit(1);
      }
      else mpilentot=(int)lentot;
      MPI_Datatype mpidtype=dtype_mpi[dtype];
      if (use_helper_server==0) {
        mpierr=mpiga_create(mpilentot, mpidtype, handle );
      }
      else {
        if (storetype==0)
         mpierr=mpiga_create(mpilentot, mpidtype, handle );
        else {
         int mproc=0;
         mpierr=twosided_helpga_create(mproc, mpilentot, handle, name, dtype);
        }
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Create: array %s created, dtype=%d, storetype=%d\n",ProcID(),name,storetype,dtype);

      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Destroy(int handle) {
      int mpierr;

      if ( mpiga_inquire_storetype(handle) == 0 )
         mpierr=mpiga_free(handle);
      else {
         int mproc=-NProcs_Work();
         mpierr=twosided_helpga_col(mproc, handle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Destroy: array %d destroyed!\n",ProcID(),handle);
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Distrib(int handle,int rank,int64_t *ilo,int64_t *ihi) {
      int mpiilo;
      int mpiihi;
      int mpierr;

      if ( mpiga_inquire_storetype(handle) == 0 ) {
         mpierr=mpiga_distribution(handle, rank, &mpiilo, &mpiihi);
      }
      else {
         mpierr=twosided_helpga_distrib(handle, rank, &mpiilo, &mpiihi);
      }

      *ilo = (int64_t) mpiilo;
      *ihi = (int64_t) mpiihi;
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Location(int handle, int64_t ilo, int64_t ihi, int64_t *map, int *proclist, int *np) {
      int mpiilo=(int)ilo;
      int mpiihi=(int)ihi;
      int mpisize;
      int mpierr;

      MPI_Comm_size(mpiga_compute_comm(), &mpisize);
      std::vector<int> mpimap(2*mpisize);
      std::vector<int> mpiproclist(mpisize);

      if ( mpiga_inquire_storetype(handle) == 0 ) {
         mpierr=mpiga_location(handle, mpiilo, mpiihi, mpimap.data(), mpiproclist.data(), np);
      }
      else {
         mpierr=twosided_helpga_location(handle, mpiilo, mpiihi, mpimap.data(), mpiproclist.data(), np);
      }

      for (int i=0;i<*np;i++) {
         map[2*i]=(int64_t)mpimap[2*i];
         map[2*i+1]=(int64_t)mpimap[2*i+1];
         proclist[i]=mpiproclist[i];
      }
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Get(int handle,int64_t ilo,int64_t ihi,void *buff) {
      int mpiilo=(int)ilo;
      int mpiihi=(int)ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(handle) == 0 )
         mpierr=mpiga_get(handle, mpiilo, mpiihi, buff);
      else {
         int mproc=0;
         int ielem=mpiilo;
         int64_t val;
         MPI_Datatype dtype=twosided_helpga_inquire_dtype(handle);
         if ( (mpiilo==mpiihi) && (dtype==MPI_INT32_T||dtype==MPI_INT64_T) ) { /* PPIDD_helpga_get_inum */
            val=twosided_helpga_one(mproc, 1, ielem, handle);
            if (dtype==MPI_INT32_T) {
              int32_t *ibuff;
              ibuff=(int32_t *)buff;
              *ibuff=(int32_t)val;
	    }
	    else {
              int64_t *ibuff;
              ibuff=(int64_t *)buff;
              *ibuff=(int64_t)val;
	    }
         }
         else if (mpiilo <= mpiihi) { /* PPIDD_helpga_get */
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra(mproc, nelem, ielem, handle,buff);
         }
         else {
            MPIGA_Error("PPIDD_Get: starting index > ending index, handle=",handle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Get: Get value from array handle= %d [%d--%d].\n",ProcID(),handle,mpiilo,mpiihi);
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Put(int handle,int64_t ilo,int64_t ihi,void *buff) {
      int mpiilo=(int)ilo;
      int mpiihi=(int)ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(handle) == 0 )
         mpierr=mpiga_put(handle, mpiilo, mpiihi, buff);
      else {
         int mproc=-NProcs_Work();
         int ielem=mpiilo;
         MPI_Datatype dtype=twosided_helpga_inquire_dtype(handle);
         if ( (mpiilo==mpiihi) && (dtype==MPI_INT32_T||dtype==MPI_INT64_T) ) { /* PPIDD_helpga_put_inum */
            int64_t valput;
            if (dtype==MPI_INT32_T) {
              int32_t *ibuff=(int32_t *)buff;
              valput=(int64_t)*ibuff;
	    }
	    else {
              int64_t *ibuff=(int64_t *)buff;
              valput=(int64_t)*ibuff;
	    }
            twosided_helpga_one(mproc, valput, ielem, handle);
         }
         else if (mpiilo <= mpiihi) { /* PPIDD_helpga_put */
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra(mproc, nelem, ielem, handle,buff);
         }
         else {
            MPIGA_Error("PPIDD_Put: starting index > ending index, handle=",handle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Put: Put buff numbers to array handle=%d [%d--%d].\n",ProcID(),handle,mpiilo,mpiihi);
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Acc(int handle,int64_t ilo,int64_t ihi,void *buff,void *fac) {
      int mpiilo=(int)ilo;
      int mpiihi=(int)ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(handle) == 0 )
         mpierr=mpiga_acc(handle, mpiilo, mpiihi, buff, fac);
      else {
         if (mpiilo <= mpiihi) { /* PPIDD_helpga_acc */
            int mproc=NProcs_Work();
            int ielem=mpiilo;
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra_acc(mproc, nelem, ielem, handle, buff, fac);
         }
         else {
            MPIGA_Error("PPIDD_Acc: starting index > ending index, handle=",handle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Acc: Accumulate buff numbers to array handle=%d [%d--%d].\n",ProcID(),handle,mpiilo,mpiihi);
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int64_t PPIDD_Read_inc(int handle,int64_t inum,int64_t incr) {
      int mpiinum = (int)inum;
      int mpiincr = (int)incr;
      int64_t mpivalue;

      if ( mpiga_inquire_storetype(handle) == 0 )
         mpivalue=(int64_t)mpiga_read_inc(handle,mpiinum,mpiincr);
      else {                                              /* PPIDD_helpga_readinc */
         int mproc=NProcs_Work();
         mpivalue=twosided_helpga_one(mproc, incr, mpiinum, handle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Read_inc: fetch-and-add element[%d] of array handle=%d by increment=%d\n",
                            ProcID(),mpiinum,handle,mpiincr);
      return (int64_t)mpivalue;
   }


   void PPIDD_Zero_patch(int handle,int64_t ilo,int64_t ihi) {
      int mpiilo = (int) ilo;
      int mpiihi = (int) ihi;
      int mpierr=0;
      if ( mpiga_inquire_storetype(handle) == 0 )
         mpierr=mpiga_zero_patch(handle,mpiilo,mpiihi);
      else
         MPIGA_Error("PPIDD_Zero_patch: invalid storetype, should be 0. handle=",handle);

      if(mpierr!=0) MPI_Abort(mpiga_compute_comm(),911);
   }


   int PPIDD_Zero(int handle) {
      int mpierr;

      if ( mpiga_inquire_storetype(handle) == 0 )
         mpierr=mpiga_zero(handle);
      else {
         int mproc=NProcs_Work();
         mpierr=twosided_helpga_col(mproc, handle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Zero: array %d has been set to zero.\n",ProcID(),handle);
      if (mpierr==0) return 1 ;
      else return 0 ;
   }


   int64_t PPIDD_Nxtval(int numproc) {
      if (use_helper_server==0) {
        fprintf(stderr,"%5d: ERROR: Attemp to call NXTVAL routine without helper process!\n", ProcID());
        MPI_Abort(mpiga_compute_comm(),911);
      }
      return (int64_t) NXTVAL(numproc);
   }


   int PPIDD_Inquire_stype(int handle) {
      return mpiga_inquire_storetype(handle);
   }


   size_t PPIDD_Inquire_mem() {
      return mpiga_localmem();
   }


   int PPIDD_Create_mutexes(int storetype,int number) {
      int mpierr;

      if (use_helper_server==0) {
        mpierr=mpiga_create_mutexes(number);      /* mutexes data store by a global array across the distributed processes */
      }
      else {
        if (storetype==0)
         mpierr=mpiga_create_mutexes(number);      /* mutexes data store by a global array across the distributed processes */
        else
         mpierr=alloc_general_helpmutexes(number);  /* mutexes data on helper process */
      }

      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   void PPIDD_Lock_mutex(int inum) {
      int mpierr;
      if ( mpigv::nmutex > 0 )
         mpierr=mpiga_lock_mutex(inum);   /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=lock_general_helpmutex(inum); /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpiga_compute_comm(),911);
   }


   void PPIDD_Unlock_mutex(int inum) {
      int mpierr;
      if ( mpigv::nmutex > 0 )
         mpierr=mpiga_unlock_mutex(inum);    /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=unlock_general_helpmutex(inum);  /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpiga_compute_comm(),911);
   }


   int PPIDD_Destroy_mutexes() {
      int mpierr;
      if ( mpigv::nmutex > 0 )
         mpierr=mpiga_destroy_mutexes();  /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=free_general_helpmutexes();   /* mutexes data on helper process */
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Eaf_open(const char *name,int type, int *handle) {
      MPI_File mpi_fh;
      MPI_Fint mpifhandle;
      MPI_Comm mpicomm=MPIGA_WORK_COMM;
      int amode=0;

      if(MPI_Debug)printf("In PPIDD_Eaf_open: begin.\n");
      switch(type){
        case MPI_EAF_RW: amode = MPI_MODE_RDWR|MPI_MODE_CREATE ;
                         break;
        case MPI_EAF_W:  amode = MPI_MODE_WRONLY|MPI_MODE_CREATE ;
                         break;
        case MPI_EAF_R:  amode = MPI_MODE_RDONLY ;
                         break;
        default:
                         MPI_Abort(mpicomm,911);
      }
/*      MPI_MODE_RDWR|MPI_MODE_CREATE|MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_UNIQUE_OPEN;*/

      int ierr=MPI_File_open(MPI_COMM_SELF,name,amode,MPI_INFO_NULL,&mpi_fh);
      mpifhandle = MPI_File_c2f( mpi_fh );
      *handle=(int)mpifhandle;
      if(MPI_Debug)printf("In PPIDD_Eaf_open: end. handle=%d,ierr=%d\n",*handle,ierr);
      return ierr;
   }


   int PPIDD_Eaf_write(int handle,double byte_offset,void *buff,size_t byte_length) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,byte_offset,(long)byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(byte_offset);
      count=(int)(byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : before MPI_File_write_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_write_at(mpi_fh,offset,buff,count,datatype,&status);
      if(MPI_Debug)printf("In PPIDD_Eaf_write : end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
   }


   int PPIDD_Eaf_awrite(int handle,double byte_offset,void *buff,size_t byte_length,int *request_id) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
      MPI_Request request;
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,byte_offset,(long)byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(byte_offset);
      count=(int)(byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : before MPI_File_iwrite_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_iwrite_at(mpi_fh,offset,buff,count,datatype,&request);
      *request_id = MPI_Request_c2f(request);
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,ierr,*request_id,(long)request);
      return ierr;
   }


   int PPIDD_Eaf_read(int handle,double byte_offset,void *buff,size_t byte_length) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,byte_offset,(long)byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(byte_offset);
      count=(int)(byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : before MPI_File_read_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_read_at(mpi_fh,offset,buff,count,datatype,&status);
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
   }


   int PPIDD_Eaf_aread(int handle,double byte_offset,void *buff,size_t byte_length,int *request_id) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
      MPI_Request request;
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,byte_offset,(long)byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(byte_offset);
      count=(int)(byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : before MPI_File_iread_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_iread_at(mpi_fh,offset,buff,count,datatype,&request);
      *request_id = MPI_Request_c2f(request);
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,ierr,*request_id,(long)request);
      return ierr;
   }


   int PPIDD_Eaf_wait(int handle,int request_id) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_Request request = MPI_Request_f2c(request_id);
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_wait  : begin. handle=%d,request_id=%d,request=%ld\n",(int)mpifhandle,request_id,(long)request);

      int ierr=MPI_Wait( &request, &status );
      if(MPI_Debug)printf("In PPIDD_Eaf_wait  : end. ierr=%d\n",ierr);
      return ierr;
   }


   int PPIDD_Eaf_waitall(const int *list, int num) {
      std::vector<MPI_Request> array_of_requests(num);
      std::vector<MPI_Status> array_of_statuses(num);
      for(int i=0;i<num;i++) array_of_requests[i] = MPI_Request_f2c(list[i]);
      return MPI_Waitall(num,array_of_requests.data(),array_of_statuses.data());
   }


   int PPIDD_Eaf_probe(int request_id,int *status) {
      int flag;
      MPI_Request request = MPI_Request_f2c(request_id);
      MPI_Status mpistatus;
      if(MPI_Debug)printf("In PPIDD_Eaf_probe  : begin. request_id=%d,request=%ld\n",request_id,(long)request);

      int ierr=MPI_Test(&request, &flag, &mpistatus);
      if(flag) *status=(int)0;
      else *status=(int)1;

      if(MPI_Debug)printf("In PPIDD_Eaf_probe  : end. ierr=%d\n",ierr);
      return ierr;
   }


   int PPIDD_Eaf_close(int handle) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      if(MPI_Debug)printf("In PPIDD_Eaf_close: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_close( &mpi_fh );
      if(MPI_Debug)printf("In PPIDD_Eaf_close: end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
   }


   int PPIDD_Eaf_delete(const char *name) {
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: begin. name=%s\n",name);
      int ierr=MPI_File_delete(name,MPI_INFO_NULL);
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: end. ierr=%d\n",ierr);
      return ierr;
   }


   int PPIDD_Eaf_length(int handle,double *fsize) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset size;
      if(MPI_Debug)printf("In PPIDD_Eaf_length: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_get_size(mpi_fh, &size);
      *fsize=(double)size;
      if(MPI_Debug)printf("In PPIDD_Eaf_length: end. handle=%d,fsize=%f\n",(int)mpifhandle,*fsize);
      return ierr;
   }


   int PPIDD_Eaf_truncate(int handle,double offset) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset size=(MPI_Offset)(offset);
      if(MPI_Debug)printf("In PPIDD_Eaf_truncate: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_set_size(mpi_fh,size);
      if(MPI_Debug)printf("In PPIDD_Eaf_truncate: end. handle=%d,size=%ld\n",(int)mpifhandle,(long)size);
      return ierr;
   }


   void PPIDD_Eaf_errmsg(int code,char *message) {
      int eclass, len;
      char estring[MPI_MAX_ERROR_STRING],estring2[MPI_MAX_ERROR_STRING];

      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: code=%d\n",code);
      MPI_Error_class(code, &eclass);
      MPI_Error_string(code, estring, &len);
      sprintf(estring2," Error %d: %s", eclass, estring);
      strncpy(message,estring2,strlen(message));
      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: message=%s\n",message);
   }


   int PPIDD_Sf_create(const char *name, double size_hard_limit, double size_soft_limit, double req_size, int *handle) {
      MPI_Comm mpicomm=MPIGA_WORK_COMM;
      MPI_File mpi_fh;

      if(MPI_Debug)printf("In PPIDD_Sf_create: begin.\n");
      int ierr=MPI_File_open(mpicomm,name,MPI_MODE_RDWR|MPI_MODE_CREATE|MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_UNIQUE_OPEN,MPI_INFO_NULL,&mpi_fh);
      MPI_Fint mpifhandle = MPI_File_c2f( mpi_fh );
      *handle=(int)mpifhandle;

      if(MPI_Debug)printf("In PPIDD_Sf_create: end. handle=%d,ierr=%d\n",*handle,ierr);
      return ierr;
   }


   int PPIDD_Sf_write(int handle, double byte_offset, double byte_length, double *buff, int *request_id) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Request request;
      if(MPI_Debug)printf("In PPIDD_Sf_write : begin. handle=%d,byte_offset=%f,byte_length=%f\n",(int)mpifhandle,byte_offset,byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(byte_offset);
      count=(int)((byte_length)/sizeof(double));
      if(MPI_Debug)printf("In PPIDD_Sf_write : before MPI_File_iwrite_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_iwrite_at(mpi_fh,offset,buff,count,MPI_DOUBLE,&request);
      *request_id = MPI_Request_c2f(request);
      if(MPI_Debug)printf("In PPIDD_Sf_write : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,ierr,*request_id,(long)request);
      return ierr;
   }


   int PPIDD_Sf_read(int handle, double byte_offset, double byte_length, double *buff, int *request_id) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Request request;
      if(MPI_Debug) printf("In PPIDD_Sf_read  : begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(byte_offset);
      count=(int)((byte_length)/sizeof(double));
      if(MPI_Debug)printf("In PPIDD_Sf_read  : before MPI_File_iread_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_iread_at(mpi_fh,offset,buff,count,MPI_DOUBLE,&request);
      *request_id = MPI_Request_c2f(request);
      if(MPI_Debug)printf("In PPIDD_Sf_read  : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,ierr,*request_id,(long)request);
      return ierr;
   }


   int PPIDD_Sf_wait(int request_id) {
      MPI_Request request = MPI_Request_f2c(request_id);
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Sf_wait  : begin. request_id=%d,request=%ld\n",request_id,(long)request);
      int ierr=MPI_Wait( &request, &status );
      if(MPI_Debug)printf("In PPIDD_Sf_wait  : end. ierr=%d\n",ierr);
      return ierr;
   }


   int PPIDD_Sf_waitall(const int *list, int num) {
      std::vector<MPI_Request> array_of_requests(num);
      std::vector<MPI_Status> array_of_statuses(num);
      for(int i=0;i<num;i++) array_of_requests[i] = MPI_Request_f2c(list[i]);
      return MPI_Waitall(num,array_of_requests.data(),array_of_statuses.data());
   }


   int PPIDD_Sf_destroy(int handle) {
      MPI_Fint mpifhandle=(MPI_Fint)handle;
      MPI_File mpi_fh;
      if(MPI_Debug)printf("In PPIDD_Sf_destroy: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_close( &mpi_fh );
      if(MPI_Debug)printf("In PPIDD_Sf_destroy: end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
   }


   void PPIDD_Sf_errmsg(int code,char *message) {
      int eclass, len;
      char estring[MPI_MAX_ERROR_STRING],estring2[MPI_MAX_ERROR_STRING];

      if(MPI_Debug)printf("In PPIDD_Sf_errmsg: code=%d\n",code);
      MPI_Error_class(code, &eclass);
      MPI_Error_string(code, estring, &len);
      sprintf(estring2," Error %d: %s", eclass, estring);
      strncpy(message,estring2,strlen(message));
      if(MPI_Debug)printf("In PPIDD_Sf_errmsg: message=%s\n",message);
   }

}
}
#endif
