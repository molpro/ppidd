#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif
#ifdef HAVE_GA_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <climits>
#include <vector>


#include <ga.h>
#include <ga-mpi.h>
#include <macdecls.h>
extern "C" {
#include <eaf.h>
#include <sf.h>
}
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

#include "mpi_utils.h"
#include "ppidd_ga_mpi.h"

static int dtype_ga(int dtype) {
 std::string errmsg;
 switch (dtype) {
  case PPIDD_FORTINT :
   if (sizeof(FORTINT)==sizeof(int)) return MT_C_INT;
   if (sizeof(FORTINT)==sizeof(long)) return MT_C_LONGINT;
   if (sizeof(FORTINT)==sizeof(long long)) return MT_C_LONGLONG;
   errmsg=" dtype_ga: unable to map FORTINT ";
   GA_Error(&errmsg[0],dtype);
   break;
  case PPIDD_DOUBLE :
   return MT_C_DBL;
  default:
   errmsg=" dtype_ga: wrong data type ";
   GA_Error(&errmsg[0],dtype);
 }
 return -1;
}

namespace ga_mpi {

   static int MPIGA_Debug=0;
   static int MPI_Debug=0;

   void PPIDD_Initialize(int *argc, char ***argv, int impl) {
    GA_Initialize_args(argc,argv);            /* initialize GA */
   }


   void PPIDD_Initialize_data() {
   }


   int64_t PPIDD_Worker_comm() {
      MPI_Comm mpicomm = GA_MPI_Comm();
      MPI_Fint fcomm=MPI_Comm_c2f(mpicomm);
      return (int64_t)fcomm;
   }


   void PPIDD_Finalize() {
      GA_Terminate();
      MPI_Finalize();
   }


   int PPIDD_Uses_ma() {
    if (GA_Uses_ma()) return 1;
    else return 0;
   }


   int PPIDD_MA_init(int dtype, int64_t stack, int64_t heap) {
      Integer istack=(Integer)stack;
      Integer iheap=(Integer)heap;
      Integer gadtype=(Integer)dtype_ga(dtype);
      if( MA_init(gadtype, istack, iheap)) return 1;
      else return 0;
   }


   double PPIDD_Wtime() {
      return GA_Wtime();
   }


   void PPIDD_Error(char *message, int code) {
      GA_Error(message,code);
   }


   void PPIDD_Helper_server(int flag, int numprocs_per_server) {
   }


   int PPIDD_Size_all() {
      int mpinp;
      MPI_Comm_size(MPI_COMM_WORLD,&mpinp);
      return mpinp;
   }


   int PPIDD_Size() {
      return GA_Nnodes();
   }


   int PPIDD_Rank() {
      return GA_Nodeid();
   }


   void PPIDD_Init_fence() {
      GA_Init_fence();
   }


   void PPIDD_Fence() {
      GA_Fence();
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
      MPI_Comm mpicomm = GA_MPI_Comm();
      int mpitag=dtype;
      int mpierr;

      int mpilenbuf = count * dtype_size(dtype);

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


   void PPIDD_Recv(void *buf,int count,int dtype,int source,int64_t *lenreal,int64_t *sourcereal,int sync) {
      MPI_Comm mpicomm = GA_MPI_Comm();
      int mpitag=dtype;
      int mpinode,mpierr;
      MPI_Status status;
      MPI_Request request;

      int mpilenbuf = count * dtype_size(dtype);

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

         *sourcereal = (int64_t) mpinode;          /* Get source node  */
         *lenreal = (int64_t) (-1);
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
         *sourcereal = (int64_t)status.MPI_SOURCE;
         *lenreal    = (int64_t)mpilenbuf;
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
      MPI_Comm mpicomm = GA_MPI_Comm();
      int flag;
      MPI_Status status;

      int mpisource = (source < 0) ? MPI_ANY_SOURCE : source;
      int mpierr = MPI_Iprobe(mpisource, tag, mpicomm, &flag, &status);
      mpi_test_status("PPIDD_Iprobe:",mpierr);
      return flag ? 1 : 0;
   }


   void PPIDD_BCast(void *buffer,int count,int dtype,int root) {
      char *cbuf=(char *)buffer;
      for (int64_t remaining=count, addr=0; remaining > 0; remaining-=(int64_t)BCAST_BATCH_SIZE, addr+=(int64_t)BCAST_BATCH_SIZE)
       GA_Brdcst(&cbuf[addr*dtype_size(dtype)], (int)(std::min((int64_t)BCAST_BATCH_SIZE,remaining)*dtype_size(dtype)) , root);
   }


   void PPIDD_Barrier() {
      GA_Sync();
   }


   void PPIDD_Gsum(int dtype,void *buffer,int len, char *op) {
      int gadtype=dtype_ga(dtype);
      GA_Gop(gadtype, buffer, len, op);
   }


   int PPIDD_Create_irreg(char *name, int64_t *lenin, int64_t nchunk, int dtype, int storetype, int *handle) {
      int ndim=1;
      ga_int nblock=(ga_int)nchunk;
      int np;
      int i;
      ga_int iad;
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
      *handle=NGA_CREATE_IRREG(gadtype, ndim, dims, name, block, &map[0]);

      return 1 ;
   }


   int PPIDD_Create(char *name,int64_t lentot, int dtype, int storetype, int *handle) {
      ga_int galentot=(ga_int)lentot;
      int gadtype=dtype_ga(dtype);

      ga_int dims[1]={galentot};
      ga_int block[1]={-1};

/*      printf("\n NGA_CREATE: %s created, dims=%d, ndim=%d\n",name,*dims,ndim); */
      *handle=NGA_CREATE(gadtype, 1, dims, name, block);

      return 1 ;
   }


   int PPIDD_Destroy(int handle) {
      GA_Destroy(handle);
      return 1 ;
   }


   int PPIDD_Distrib(int handle,int rank,int64_t *ilo,int64_t *ihi) {
      ga_int gailo[1];
      ga_int gaihi[1];

      NGA_DISTRIBUTION(handle, rank, gailo, gaihi);
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
   }


   int PPIDD_Location(int handle, int64_t *ilo, int64_t *ihi, int64_t *map, int64_t *proclist, int *np) {
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      int mpisize = GA_Nnodes();
      std::vector<ga_int> mpimap(2*mpisize);
      std::vector<int> mpiproclist(mpisize);

      *np=NGA_LOCATE_REGION(handle, mpiilo, mpiihi, &mpimap[0], &mpiproclist[0]);

      for (int i=0;i<*np;i++) {
         map[2*i]=(int64_t)(mpimap[2*i]+1);
         map[2*i+1]=(int64_t)(mpimap[2*i+1]+1);
         proclist[i]=(int64_t)mpiproclist[i];
      }
      return 1 ;
   }


   int PPIDD_Get(int handle,int64_t *ilo,int64_t *ihi,void *buff) {
      ga_int ld[1]={1};
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_GET(handle, mpiilo, mpiihi, buff, ld);
      return 1 ;
   }


   int PPIDD_Put(int handle,int64_t *ilo,int64_t *ihi,void *buff) {
      ga_int ld[1]={1};
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_PUT(handle, mpiilo, mpiihi, buff, ld);
      return 1 ;
   }


   int PPIDD_Acc(int handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac) {
      ga_int ld[1]={1};
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_ACC(handle, mpiilo, mpiihi, buff, ld, fac);
      return 1 ;
   }


   int64_t PPIDD_Read_inc(int handle,int64_t inum,int64_t incr) {
      ga_int mpiinum[1];

      mpiinum[0] = (ga_int) inum-1;
      return (int64_t) NGA_READ_INC(handle,mpiinum, (long)incr);
   }


   void PPIDD_Zero_patch(int handle,int64_t *ilo,int64_t *ihi) {
      ga_int mpiilo[1]={(ga_int)*ilo-1};
      ga_int mpiihi[1]={(ga_int)*ihi-1};

      NGA_ZERO_PATCH(handle, mpiilo, mpiihi);
   }


   int PPIDD_Zero(int handle) {
      GA_Zero(handle);
      return 1 ;
   }


   static int PPIDD_Nxtval_initialised=0;
   static int PPIDD_Nxtval_handle;
   void PPIDD_Nxtval(int numproc, int64_t *val) {
      if (numproc < 0) {
        /* reset - collective */
        if (PPIDD_Nxtval_initialised) PPIDD_Destroy(PPIDD_Nxtval_handle);
        PPIDD_Nxtval_initialised=0;
        //      }
        //else if (! PPIDD_Nxtval_initialised) {
        /* first call needs to be collective and will return 0*/
	std::string name="Nxtval";
        PPIDD_Create(&name[0],1,0,1,&PPIDD_Nxtval_handle);
        PPIDD_Zero(PPIDD_Nxtval_handle);
        PPIDD_Nxtval_initialised=1;
        *val=0;
      }
      else {
        *val = PPIDD_Read_inc(PPIDD_Nxtval_handle,1,1);
      }
   }


   void PPIDD_Inquire_name(int handle, char *name) {
      strncpy(name,GA_Inquire_name(handle),strlen(name));
   }


   int PPIDD_Inquire_stype(int handle) {
      return 0;
   }


   void PPIDD_Inquire_mem(int64_t *mem_used) {
      *mem_used=(int64_t)GA_Inquire_memory();
   }


   int PPIDD_Create_mutexes(int storetype,int number) {
      int mpierr = GA_Create_mutexes(number);
/* This is one of exceptions in GA (see global/src/capi.c) : Returns [1] if the operation succeeded or [0] when failed */
      if(MPIGA_Debug)printf("In PPIDD_Create_Mutexes: mpierr=%d.\n",mpierr);
      if(mpierr==1) return 1 ;
      else return 0 ;
   }


   void PPIDD_Lock_mutex(int inum) {
      GA_Lock(inum);
   }


   void PPIDD_Unlock_mutex(int inum) {
      GA_Unlock(inum);
   }


   int PPIDD_Destroy_mutexes() {
      int mpierr=GA_Destroy_mutexes();
/* This is one of exceptions in GA (see global/src/capi.c) : Returns [1] if the operation succeeded or [0] when failed */
      if(MPIGA_Debug)printf("In PPIDD_Destroy_Mutexes: mpierr=%d.\n",mpierr);
      if(mpierr==1) return 1 ;
      else return 0 ;
   }


   int PPIDD_Eaf_open(char *name,int type, int *handle) {
      if(MPI_Debug)printf("In PPIDD_Eaf_open: begin.\n");
      int ierr=EAF_Open(name, type, handle);
      if(MPI_Debug)printf("In PPIDD_Eaf_open: end. handle=%d,ierr=%d\n",*handle,ierr);
      return ierr;
   }


   int PPIDD_Eaf_write(int handle,double *byte_offset,void *buff,int64_t *byte_length) {
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;

      return EAF_Write(handle,offset,buff,bytes);
   }


   int PPIDD_Eaf_awrite(int handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int request;

      int ierr=EAF_Awrite(handle,offset,buff,bytes,&request);
      *request_id=(int64_t)request;
      return ierr;
   }


   int PPIDD_Eaf_read(int handle,double *byte_offset,void *buff,int64_t *byte_length) {
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;

      return EAF_Read(handle,offset,buff,bytes);
   }


   int PPIDD_Eaf_aread(int handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int request;

      int ierr=EAF_Aread(handle,offset,buff,bytes,&request);
      *request_id=(int64_t)request;
      return ierr;
   }


   int PPIDD_Eaf_wait(int handle,int64_t *request_id) {
      int request=(int)*request_id;

      return EAF_Wait(handle,request);
   }


   int PPIDD_Eaf_waitall(int64_t *list, int num) {
      return 0;
   }


   int PPIDD_Eaf_probe(int64_t *request_id,int *status) {
      int garequest=(int)*request_id;

      return EAF_Probe(garequest, status);
   }


   int PPIDD_Eaf_close(int handle) {
      return EAF_Close(handle);
   }


   int PPIDD_Eaf_delete(char *name) {
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: begin. name=%s\n",name);
      int ierr=EAF_Delete(name);
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: end. ierr=%d\n",ierr);
      return ierr;
   }


   int PPIDD_Eaf_length(int handle,double *fsize) {
      eaf_off_t length;

      int ierr=EAF_Length(handle,&length);
      *fsize=(double)length;
      return ierr;
   }


   int PPIDD_Eaf_truncate(int handle,double *offset) {
      eaf_off_t length=(eaf_off_t)*offset;

      return EAF_Truncate(handle,length);
   }


   void PPIDD_Eaf_errmsg(int code,char *message) {
      EAF_Errmsg(code, message);
   }


   int PPIDD_Sf_create(char *name ,double *size_hard_limit, double *size_soft_limit, double *req_size, int *handle) {

      if(MPI_Debug)printf("In PPIDD_Sf_create: begin.\n");
      int ierr=SF_Create(name, *size_hard_limit, *size_soft_limit, *req_size, handle);

      if(MPI_Debug)printf("In PPIDD_Sf_create: end. handle=%d,ierr=%d\n",*handle,ierr);
      return ierr;
   }


   int PPIDD_Sf_write(int handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
      char *buffer=(char *)buff;

      int irequest_id;
      int ierr=SF_Write(handle, *byte_offset, *byte_length, buffer, &irequest_id);
      *request_id=(int64_t)irequest_id;
      return ierr;
   }


   int PPIDD_Sf_read(int handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
      char *buffer=(char *)buff;

      int irequest_id;
      int ierr=SF_Read(handle, *byte_offset, *byte_length, buffer, &irequest_id);
      *request_id=(int64_t)irequest_id;
      return ierr;
   }


   int PPIDD_Sf_wait(int64_t *request_id) {
      int irequest_id=(int)*request_id;
      int ierr=SF_Wait(&irequest_id);
      *request_id=(int64_t)irequest_id;
      return ierr;
   }


   int PPIDD_Sf_waitall(int64_t *list, int num) {
      int *ilist;

      ilist=(int *)malloc(num * sizeof(int));
      for (int i=0; i<num; ++i) ilist[i] = (int) list[i];

      int ierr=SF_Waitall(ilist, num);

      for (int i=0; i<num; ++i) list[i] = (int64_t) ilist[i];
      free(ilist);
      return ierr;
   }


   int PPIDD_Sf_destroy(int handle) {
      return SF_Destroy(handle);
   }


   void PPIDD_Sf_errmsg(int code,char *message) {
      SF_Errmsg(code, message);
   }

}
#endif
