#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include "ppidd_no_mpi.h"

namespace no_mpi {

   [[ noreturn ]] static void do_not_call(const char* function) {
    fprintf(stderr,"%s should not be called in no_mpi case\n",function);
    exit(1);
   }

   void PPIDD_Initialize(int *argc, char ***argv, int impl) {
   }


   void PPIDD_Initialize_data() {
   }


   int64_t PPIDD_Worker_comm() {
      do_not_call("PPIDD_Worker_comm");
   }


   void PPIDD_Finalize() {
   }


   int PPIDD_Uses_ma() {
      return 0;
   }


   int PPIDD_MA_init(int dtype, int64_t *stack, int64_t *heap) {
      return 1;
   }


   void PPIDD_Wtime(double *ctime) {
      *ctime = (double)0;
   }


   void PPIDD_Error(char *message, int code) {
      fprintf(stdout," %s %d (%#x).\n", message,code,code);
      fflush(stdout);
      fprintf(stderr," %s %d (%#x).\n", message,code,code);

      printf(" PPIDD_Error: now exiting...\n");
      exit(1);
   }


   void PPIDD_Helper_server(int *flag, int64_t *numprocs_per_server) {
   }


   void PPIDD_Size_all(int64_t *np) {
      *np = (int64_t)1;
   }


   void PPIDD_Size(int64_t *np) {
      *np = (int64_t)1;
   }


   void PPIDD_Rank(int64_t *me) {
      *me = (int64_t)0;
   }


   void PPIDD_Init_fence() {
   }


   void PPIDD_Fence() {
   }


   void PPIDD_Send(void *buf,int64_t *count,int dtype,int64_t *dest,int64_t *sync) {
      do_not_call("PPIDD_Send");
   }


   void PPIDD_Recv(void *buf,int64_t *count,int dtype,int64_t *source,int64_t *lenreal,int64_t *sourcereal,int64_t *sync) {
      do_not_call("PPIDD_Recv");
   }


   void PPIDD_Wait(int64_t *nodesel) {
   }


   int PPIDD_Iprobe(int64_t *tag,int64_t *source) {
    return 0 ;
   }


   void PPIDD_BCast(void *buffer,int64_t *count,int dtype,int64_t *root) {
   }


   void PPIDD_Barrier() {
   }


   void PPIDD_Gsum(int dtype,void *buffer,int64_t *len, char *op) {
   }


   int PPIDD_Create_irreg(char *name, int64_t *lenin, int64_t *nchunk, int dtype, int64_t *storetype, int64_t *handle) {
      do_not_call("PPIDD_Create_irreg");
   }


   int PPIDD_Create(char *name,int64_t *lentot, int dtype, int64_t *storetype, int64_t *handle) {
      do_not_call("PPIDD_Create");
   }


   int PPIDD_Destroy(int64_t *handle) {
      do_not_call("PPIDD_Destroy");
   }


   int PPIDD_Distrib(int64_t *handle,int64_t *rank,int64_t *ilo,int64_t *ihi) {
      do_not_call("PPIDD_Distrib");
   }


   int PPIDD_Location(int64_t *handle, int64_t *ilo, int64_t *ihi, int64_t *map, int64_t *proclist, int64_t *np) {
      do_not_call("PPIDD_Location");
   }


   int PPIDD_Get(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
      do_not_call("PPIDD_Get");
   }


   int PPIDD_Put(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff) {
      do_not_call("PPIDD_Put");
   }


   int PPIDD_Acc(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac) {
      do_not_call("PPIDD_Acc");
   }


   void PPIDD_Read_inc(int64_t *handle,int64_t *inum,int64_t *incr,int64_t *returnval) {
      do_not_call("PPIDD_Read_inc");
   }


   void PPIDD_Zero_patch(int64_t *handle,int64_t *ilo,int64_t *ihi) {
      do_not_call("PPIDD_Zero_patch");
   }


   int PPIDD_Zero(int64_t *handle) {
      do_not_call("PPIDD_Zero");
   }


   void PPIDD_Nxtval(int64_t *numproc, int64_t *val) {
      do_not_call("PPIDD_Nxtval");
   }


   void PPIDD_Duplicate(int64_t *handlei, int64_t *handlej, char *name) {
      do_not_call("PPIDD_Duplicate");
   }


   void PPIDD_Inquire_name(int64_t *handle, char *name) {
      do_not_call("PPIDD_Inquire_name");
   }


   void PPIDD_Inquire_stype(int64_t *handle, int64_t *storetype) {
      do_not_call("PPIDD_Inquire_stype");
   }


   void PPIDD_Inquire_mem(int64_t *mem_used) {
      *mem_used=(int64_t)0;
   }


   int PPIDD_Create_mutexes(int64_t *storetype,int64_t *number) {
      return 1 ;
   }


   void PPIDD_Lock_mutex(int64_t *inum) {
   }


   void PPIDD_Unlock_mutex(int64_t *inum) {
   }


   int PPIDD_Destroy_mutexes() {
      return 1 ;
   }


   int PPIDD_Eaf_open(char *name,int64_t *type, int64_t *handle) {
      do_not_call("PPIDD_Eaf_open");
   }


   int PPIDD_Eaf_write(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length) {
      do_not_call("PPIDD_Eaf_write");
   }


   int PPIDD_Eaf_awrite(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
      do_not_call("PPIDD_Eaf_awrite");
   }


   int PPIDD_Eaf_read(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length) {
      do_not_call("PPIDD_Eaf_read");
   }


   int PPIDD_Eaf_aread(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
      do_not_call("PPIDD_Eaf_aread");
   }


   int PPIDD_Eaf_wait(int64_t *handle,int64_t *request_id) {
      do_not_call("PPIDD_Eaf_wait");
   }


   int PPIDD_Eaf_waitall(int64_t *list, int64_t *num) {
      do_not_call("PPIDD_Eaf_waitall");
   }


   int PPIDD_Eaf_probe(int64_t *request_id,int64_t *status) {
      do_not_call("PPIDD_Eaf_probe");
   }


   int PPIDD_Eaf_close(int64_t *handle) {
      do_not_call("PPIDD_Eaf_close");
   }


   int PPIDD_Eaf_delete(char *name) {
      do_not_call("PPIDD_Eaf_delete");
   }


   int PPIDD_Eaf_length(int64_t *handle,double *fsize) {
      do_not_call("PPIDD_Eaf_length");
   }


   int PPIDD_Eaf_truncate(int64_t *handle,double *offset) {
      do_not_call("PPIDD_Eaf_truncate");
   }


   void PPIDD_Eaf_errmsg(int *code,char *message) {
      do_not_call("PPIDD_Eaf_errmsg");
   }


   int PPIDD_Sf_create(char *name ,double *size_hard_limit, double *size_soft_limit, double *req_size, int64_t *handle) {
      do_not_call("PPIDD_Sf_create");
   }


   int PPIDD_Sf_write(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
      do_not_call("PPIDD_Sf_write");
   }


   int PPIDD_Sf_read(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
      do_not_call("PPIDD_Sf_read");
   }


   int PPIDD_Sf_wait(int64_t *request_id) {
      do_not_call("PPIDD_Sf_wait");
   }


   int PPIDD_Sf_waitall(int64_t *list, int64_t *num) {
      do_not_call("PPIDD_Sf_waitall");
   }


   int PPIDD_Sf_destroy(int64_t *handle) {
      do_not_call("PPIDD_Sf_destroy");
   }


   void PPIDD_Sf_errmsg(int *code,char *message) {
      do_not_call("PPIDD_Sf_errmsg");
   }

}
