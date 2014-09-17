#ifndef __PPIDD_C_H__
#define __PPIDD_C_H__
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
   extern void PPIDD_Initialize(int argc, char **argv);
   extern void PPIDD_Initialize_data(void);
   extern int64_t  PPIDD_Worker_comm(void);
   extern void PPIDD_Finalize(void);
   extern void PPIDD_Uses_ma(int64_t *ok);
   extern void PPIDD_MA_init(int64_t *dtype, int64_t *stack, int64_t *heap, int64_t *ok);
   extern void PPIDD_Wtime(double *ctime);
   extern void PPIDD_Error(char *message,int64_t *code);
   extern void PPIDD_Helper_server(int64_t *flag, int64_t *num_per_server);
   extern void PPIDD_Nxtval(int64_t *numproc, int64_t *val);
   extern void PPIDD_Size_all(int64_t *np);
   extern void PPIDD_Size(int64_t *np);
   extern void PPIDD_Rank(int64_t *me);
   extern void PPIDD_Init_fence(void);
   extern void PPIDD_Fence(void);
   extern void PPIDD_Send(void *buf,int64_t *count,int64_t *dtype,int64_t *dest,int64_t *sync);
   extern void PPIDD_Recv(void *buf,int64_t *count,int64_t *dtype,int64_t *source,int64_t *lenreal,int64_t *sourcereal,int64_t *sync);
   extern void PPIDD_Wait(int64_t *nodesel);
   extern void PPIDD_Iprobe(int64_t *tag,int64_t *source,int64_t *ok);
   extern void PPIDD_BCast(void *buffer,int64_t *count,int64_t *type,int64_t *root);
   extern void PPIDD_Barrier(void);
   extern void PPIDD_Gsum(int64_t *type,void *buffer,int64_t *len, char *op);
   extern void PPIDD_Create_irreg(char *name,int64_t *lenin,int64_t *nchunk,int64_t *datatype,int64_t *storetype,int64_t *handle,int64_t *ok);
   extern void PPIDD_Create(char *name,int64_t *lentot, int64_t *datatype, int64_t *storetype, int64_t *handle, int64_t *ok);
   extern void PPIDD_Destroy(int64_t *handle,int64_t *ok);
   extern void PPIDD_Distrib(int64_t *handle,int64_t *rank,int64_t *ilo,int64_t *ihi,int64_t *ok);
   extern void PPIDD_Location(int64_t *handle,int64_t *ilo,int64_t *ihi,int64_t *map,int64_t *proclist,int64_t *np,int64_t *ok);
   extern void PPIDD_Get(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,int64_t *ok);
   extern void PPIDD_Put(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,int64_t *ok);
   extern void PPIDD_Acc(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac,int64_t *ok);
   extern void PPIDD_Read_inc(int64_t *ihandle,int64_t *inum,int64_t *incr,int64_t *returnval);
   extern void PPIDD_Zero_patch(int64_t *ihandle,int64_t *ilo,int64_t *ihi);
   extern void PPIDD_Zero(int64_t *handle,int64_t *ok);
   extern void PPIDD_Duplicate(int64_t *handlei, int64_t *handlej, char *name);
   extern void PPIDD_Inquire_name(int64_t *handle, char *name);
   extern void PPIDD_Inquire_stype(int64_t *handle, int64_t *storetype);
   extern void PPIDD_Inquire_mem(int64_t *mem_used);
   extern void PPIDD_Create_mutexes(int64_t *storetype,int64_t *number,int64_t *ok);
   extern void PPIDD_Lock_mutex(int64_t *inum);
   extern void PPIDD_Unlock_mutex(int64_t *inum);
   extern void PPIDD_Destroy_mutexes(int64_t *ok);
#ifdef __cplusplus
}
#endif
#endif
