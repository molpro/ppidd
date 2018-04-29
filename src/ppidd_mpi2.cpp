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
#include "ppidd_mpi2.h"


namespace mpi2 {

    static int MPIGA_Debug=0;

   static void do_not_call(const char* function) {
    fprintf(stderr,"%s should not be called in mpi2 case\n",function);
    exit(1);
   }

   void PPIDD_Initialize(int *argc, char ***argv) {
    int mpierr=mpiga_initialize(argc,argv);
    mpi_test_status("PPIDD_Initialize:",mpierr);
   }


   void PPIDD_Initialize_data(void) {
      mpiga_initialize_data();
   }


   int64_t PPIDD_Worker_comm(void) {
      MPI_Comm mycomm=mpiga_compute_comm();

/* test whether worker communicator contains all the processes, if so then return MPI_COMM_WORLD */
      int np_all, np_worker=mpigv(nprocs);
      MPI_Comm_size(MPI_COMM_WORLD, &np_all);
      if(np_all==np_worker) mycomm=MPI_COMM_WORLD;

      MPI_Fint fcomm=MPI_Comm_c2f(mycomm);
      return (int64_t)fcomm;
   }


   void PPIDD_Finalize(void) {
      mpiga_terminate();
   }


   int PPIDD_Uses_ma() {
    return 0;
   }


   int PPIDD_MA_init(int dtype, int64_t *stack, int64_t *heap) {
      return 1;
   }


   void PPIDD_Wtime(double *ctime) {
      *ctime = MPI_Wtime();
   }


   void PPIDD_Error(char *message,int *code) {
      MPIGA_Error(message,*code);
   }


   void PPIDD_Helper_server(int *flag, int64_t *numprocs_per_server) {
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
   }


   void PPIDD_Size_all(int64_t *np) {
      int mpinp;
      MPI_Comm_size(MPI_COMM_WORLD,&mpinp);
      *np = (int64_t) mpinp;
   }


   void PPIDD_Size(int64_t *np) {
      int mpinp;
      MPI_Comm_size(mpiga_compute_comm(),&mpinp);
      *np = (int64_t) mpinp;
   }


   void PPIDD_Rank(int64_t *me) {
      int mpime;
      MPI_Comm_rank(mpiga_compute_comm(),&mpime);
      *me = (int64_t) mpime;
   }


   void PPIDD_Init_fence(void) {
   }


   void PPIDD_Fence(void) {
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


   void PPIDD_Send(void *buf,int64_t *count,int dtype,int64_t *dest,int64_t *sync) {
      MPI_Comm mpicomm=mpiga_compute_comm();
      int mpicount=(int)*count;
      int mpidest=(int)*dest;
      int mpitag=dtype;
      int mpisync=(int)*sync;
      int mpierr;

      int mpilenbuf = mpicount * dtype_size(dtype);

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
   }


   void PPIDD_Recv(void *buf,int64_t *count,int dtype,int64_t *source,int64_t *lenreal,int64_t *sourcereal,int64_t *sync) {
      MPI_Comm mpicomm=mpiga_compute_comm();
      int mpicount=(int)*count;
      int mpitag=dtype;
      int mpisource=(int)*source;
      int mpisync=(int)*sync;
      int mpinode,mpierr;
      MPI_Status status;
      MPI_Request request;

      int mpilenbuf = mpicount * dtype_size(dtype);

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
   }


   void PPIDD_Wait(int64_t *nodesel) {
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


   int PPIDD_Iprobe(int64_t *tag,int64_t *source) {
      MPI_Comm mpicomm=mpiga_compute_comm();
      int mpitag=(int)*tag;
      int flag;
      MPI_Status status;

      int mpisource = (*source < 0) ? MPI_ANY_SOURCE  : (int) *source;
      int mpierr = MPI_Iprobe(mpisource, mpitag, mpicomm, &flag, &status);
      mpi_test_status("PPIDD_Iprobe:",mpierr);
      if(flag) return 1 ;
      else return 0 ;
   }


   void PPIDD_BCast(void *buffer,int64_t *count,int dtype,int64_t *root) {
      int mpicount=(int)*count;
      int mpiroot=(int)*root;

      MPI_Datatype mpidtype=dtype_mpi(dtype);

      int mpierr=MPI_Bcast(buffer,mpicount,mpidtype,mpiroot,mpiga_compute_comm());
      mpi_test_status("PPIDD_BCast:",mpierr);
   }


   void PPIDD_Barrier(void) {
      int mpierr=MPI_Barrier(mpiga_compute_comm());
      mpi_test_status("PPIDD_Barrier:",mpierr);
   }


   void PPIDD_Gsum(int dtype,void *buffer,int64_t *len, char *op) {
      int mpilen=(int)*len;
      MPI_Datatype mpidtype=dtype_mpi(dtype);
      MPI_GSum(mpidtype,buffer,mpilen, op);
   }


   int PPIDD_Create_irreg(char *name, int64_t *lenin, int64_t *nchunk, int dtype, int64_t *storetype, int64_t *handle) {
      int mpierr;
      int mpinchunk=(int)*nchunk;
      int stype=(int)*storetype;
      int mpihandle;

      std::vector<int> mpilenin(mpinchunk);
      for (int i=0;i<mpinchunk;i++) mpilenin[i]=(int)lenin[i];
      MPI_Datatype mpidtype=dtype_mpi(dtype);
      if (use_helper_server==0) {
        mpierr=mpiga_create_irreg(name, &mpilenin[0], mpinchunk, mpidtype, &mpihandle);
      }
      else {
        if (stype==0)
          mpierr=mpiga_create_irreg(name, &mpilenin[0], mpinchunk, mpidtype, &mpihandle);
        else {
          int mproc=0;
          mpierr=twosided_helpga_create_irreg(mproc, &mpilenin[0], mpinchunk, &mpihandle, name, dtype);
        }
      }
      *handle=(int64_t)mpihandle;
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Create(char *name,int64_t *lentot, int dtype, int64_t *storetype, int64_t *handle) {
      int mpierr;
      int mpilentot;
      int stype=(int)*storetype;
      int mpihandle;

      if (*lentot > INT_MAX) {
       printf(" ERROR: PPIDD_Create: lentot too large for MPI\n");
       exit(1);
      }
      else mpilentot=(int)*lentot;
      MPI_Datatype mpidtype=dtype_mpi(dtype);
      if (use_helper_server==0) {
        mpierr=mpiga_create( name, mpilentot, mpidtype, &mpihandle );
      }
      else {
        if (stype==0)
         mpierr=mpiga_create( name, mpilentot, mpidtype, &mpihandle );
        else {
         int mproc=0;
         mpierr=twosided_helpga_create(mproc, mpilentot, &mpihandle, name, dtype);
        }
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Create: array %s created, dtype=%d, storetype=%d\n",ProcID(),name,stype,dtype);

      *handle=(int64_t)mpihandle;
      if(mpierr==0) return 1 ;
      else return 0 ;
   }


   int PPIDD_Destroy(int64_t *handle) {
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
   }


   int PPIDD_Distrib(int64_t *handle,int64_t *rank,int64_t *ilo,int64_t *ihi) {
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
   }


   int PPIDD_Location(int64_t *handle, int64_t *ilo, int64_t *ihi, int64_t *map, int64_t *proclist, int64_t *np) {
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpisize,mpinp;
      int mpierr;

      MPI_Comm_size(mpiga_compute_comm(), &mpisize);
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
   }


   int PPIDD_Get(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
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
   }


   int PPIDD_Put(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
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
   }


   int PPIDD_Acc(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac) {
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
   }


   void PPIDD_Read_inc(int64_t *handle,int64_t *inum,int64_t *incr,int64_t *returnval) {
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
   }


   void PPIDD_Zero_patch(int64_t *handle,int64_t *ilo,int64_t *ihi) {
      int mpihandle = (int) *handle;
      int mpiilo = (int) *ilo;
      int mpiihi = (int) *ihi;
      int mpierr=0;
      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_zero_patch(mpihandle,mpiilo,mpiihi);
      else
         MPIGA_Error("PPIDD_Zero_patch: invalid storetype, should be 0. handle=",mpihandle);

      if(mpierr!=0) MPI_Abort(mpiga_compute_comm(),911);
   }


   int PPIDD_Zero(int64_t *handle) {
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
   }


   void PPIDD_Nxtval(int64_t *numproc, int64_t *val) {
      if (use_helper_server==0) {
        fprintf(stderr,"%5d: ERROR: Attemp to call NXTVAL routine without helper process!\n", ProcID());
        MPI_Abort(mpiga_compute_comm(),911);
      }
      else {
        int mproc = (int) *numproc;
        *val= (int64_t) NXTVAL(&mproc);
      }
   }


   void PPIDD_Duplicate(int64_t *handlei, int64_t *handlej, char *name) {
      do_not_call("PPIDD_Duplicate");
   }


   void PPIDD_Inquire_name(int64_t *handle, char *name) {
      char *name2;
      int mpihandle = (int) *handle;

      int lxi=strlen(name);
      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpiga_inquire_name(mpihandle, &name2);
      else {
         twosided_helpga_inquire_name(mpihandle, &name2);
      }
      int len_actual=strlen(name2);
      strncpy(name,name2,len_actual);
      for(int i=len_actual;i<lxi;i++) name[i]=' ';
      if(MPIGA_Debug)printf("In PPIDD_Inquire_name: name2=%s,strlen(name2)=%d,lxi=%d\n",name2,len_actual,lxi);
   }


   void PPIDD_Inquire_stype(int64_t *handle, int64_t *storetype) {
      int mpihandle = (int) *handle;
      *storetype=(int64_t)mpiga_inquire_storetype(mpihandle);
   }


   void PPIDD_Inquire_mem(int64_t *mem_used) {
      long localmem;
      localmem=mpiga_localmem();
      *mem_used=(int64_t)localmem;
   }


   int PPIDD_Create_mutexes(int64_t *storetype,int64_t *number) {
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
   }


   void PPIDD_Lock_mutex(int64_t *inum) {
      int mpiinum = (int) *inum;
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_lock_mutex(mpiinum);   /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=lock_general_helpmutex(mpiinum); /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpiga_compute_comm(),911);
   }


   void PPIDD_Unlock_mutex(int64_t *inum) {
      int mpiinum = (int) *inum;
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_unlock_mutex(mpiinum);    /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=unlock_general_helpmutex(mpiinum);  /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpiga_compute_comm(),911);
   }


   int PPIDD_Destroy_mutexes() {
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_destroy_mutexes();  /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=free_general_helpmutexes();   /* mutexes data on helper process */
      if(mpierr==0) return 1 ;
      else return 0 ;
   }

}
#endif
