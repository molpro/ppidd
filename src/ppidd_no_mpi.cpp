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


   int PPIDD_MA_init(int dtype, int64_t stack, int64_t heap) {
      return 1;
   }


   double PPIDD_Wtime() {
      return (double)0;
   }


   void PPIDD_Error(char *message, int code) {
      fprintf(stdout," %s %d (%#x).\n", message,code,code);
      fflush(stdout);
      fprintf(stderr," %s %d (%#x).\n", message,code,code);

      printf(" PPIDD_Error: now exiting...\n");
      exit(1);
   }


   void PPIDD_Helper_server(int flag, int numprocs_per_server) {
   }


   int PPIDD_Size_all() {
      return 1;
   }


   int PPIDD_Size() {
      return 1;
   }


   int PPIDD_Rank() {
      return 0;
   }


   void PPIDD_Init_fence() {
   }


   void PPIDD_Fence() {
   }


   void PPIDD_Send(void *buf,int count,int dtype,int dest,int sync) {
      do_not_call("PPIDD_Send");
   }


   void PPIDD_Recv(void *buf,int count,int dtype,int source,int64_t *lenreal,int64_t *sourcereal,int sync) {
      do_not_call("PPIDD_Recv");
   }


   void PPIDD_Wait() {
   }


   int PPIDD_Iprobe(int tag, int source) {
    return 0 ;
   }


   void PPIDD_BCast(void *buffer,int count,int dtype,int root) {
   }


   void PPIDD_Barrier() {
   }


   void PPIDD_Gsum(int dtype,void *buffer,int len, char *op) {
   }


   int PPIDD_Create_irreg(char *name, int64_t *lenin, int64_t nchunk, int dtype, int storetype, int *handle) {
      do_not_call("PPIDD_Create_irreg");
   }


   int PPIDD_Create(char *name,int64_t lentot, int dtype, int storetype, int *handle) {
      do_not_call("PPIDD_Create");
   }


   int PPIDD_Destroy(int handle) {
      do_not_call("PPIDD_Destroy");
   }


   int PPIDD_Distrib(int handle,int rank,int64_t *ilo,int64_t *ihi) {
      do_not_call("PPIDD_Distrib");
   }


   int PPIDD_Location(int handle, int64_t *ilo, int64_t *ihi, int64_t *map, int64_t *proclist, int *np) {
      do_not_call("PPIDD_Location");
   }


   int PPIDD_Get(int handle,int64_t *ilo,int64_t *ihi,void *buff) {
      do_not_call("PPIDD_Get");
   }


   int PPIDD_Put(int handle,int64_t *ilo,int64_t *ihi,void *buff) {
      do_not_call("PPIDD_Put");
   }


   int PPIDD_Acc(int handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac) {
      do_not_call("PPIDD_Acc");
   }


   int64_t PPIDD_Read_inc(int handle,int64_t inum,int64_t incr) {
      do_not_call("PPIDD_Read_inc");
   }


   void PPIDD_Zero_patch(int handle,int64_t *ilo,int64_t *ihi) {
      do_not_call("PPIDD_Zero_patch");
   }


   int PPIDD_Zero(int handle) {
      do_not_call("PPIDD_Zero");
   }


   void PPIDD_Nxtval(int numproc, int64_t *val) {
      do_not_call("PPIDD_Nxtval");
   }


   void PPIDD_Inquire_name(int handle, char *name) {
      do_not_call("PPIDD_Inquire_name");
   }


   int PPIDD_Inquire_stype(int handle) {
      do_not_call("PPIDD_Inquire_stype");
   }


   void PPIDD_Inquire_mem(int64_t *mem_used) {
      *mem_used=(int64_t)0;
   }


   int PPIDD_Create_mutexes(int storetype,int number) {
      return 1 ;
   }


   void PPIDD_Lock_mutex(int inum) {
   }


   void PPIDD_Unlock_mutex(int inum) {
   }


   int PPIDD_Destroy_mutexes() {
      return 1 ;
   }


   int PPIDD_Eaf_open(char *name,int type, int *handle) {
      do_not_call("PPIDD_Eaf_open");
   }


   int PPIDD_Eaf_write(int handle,double *byte_offset,void *buff,int64_t *byte_length) {
      do_not_call("PPIDD_Eaf_write");
   }


   int PPIDD_Eaf_awrite(int handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
      do_not_call("PPIDD_Eaf_awrite");
   }


   int PPIDD_Eaf_read(int handle,double *byte_offset,void *buff,int64_t *byte_length) {
      do_not_call("PPIDD_Eaf_read");
   }


   int PPIDD_Eaf_aread(int handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
      do_not_call("PPIDD_Eaf_aread");
   }


   int PPIDD_Eaf_wait(int handle,int64_t *request_id) {
      do_not_call("PPIDD_Eaf_wait");
   }


   int PPIDD_Eaf_waitall(int64_t *list, int num) {
      do_not_call("PPIDD_Eaf_waitall");
   }


   int PPIDD_Eaf_probe(int64_t *request_id,int* status) {
      do_not_call("PPIDD_Eaf_probe");
   }


   int PPIDD_Eaf_close(int handle) {
      do_not_call("PPIDD_Eaf_close");
   }


   int PPIDD_Eaf_delete(char *name) {
      do_not_call("PPIDD_Eaf_delete");
   }


   int PPIDD_Eaf_length(int handle,double *fsize) {
      do_not_call("PPIDD_Eaf_length");
   }


   int PPIDD_Eaf_truncate(int handle,double *offset) {
      do_not_call("PPIDD_Eaf_truncate");
   }


   void PPIDD_Eaf_errmsg(int code,char *message) {
      do_not_call("PPIDD_Eaf_errmsg");
   }


   int PPIDD_Sf_create(char *name ,double *size_hard_limit, double *size_soft_limit, double *req_size, int *handle) {
      do_not_call("PPIDD_Sf_create");
   }


   int PPIDD_Sf_write(int handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
      do_not_call("PPIDD_Sf_write");
   }


   int PPIDD_Sf_read(int handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id) {
      do_not_call("PPIDD_Sf_read");
   }


   int PPIDD_Sf_wait(int64_t *request_id) {
      do_not_call("PPIDD_Sf_wait");
   }


   int PPIDD_Sf_waitall(int64_t *list, int num) {
      do_not_call("PPIDD_Sf_waitall");
   }


   int PPIDD_Sf_destroy(int handle) {
      do_not_call("PPIDD_Sf_destroy");
   }


   void PPIDD_Sf_errmsg(int code,char *message) {
      do_not_call("PPIDD_Sf_errmsg");
   }

}
