 extern void PPIDD_Initialize(int *argc, char ***argv, int impl);
 extern void PPIDD_Initialize_data();
 extern int64_t  PPIDD_Worker_comm();
 extern void PPIDD_Finalize();
 extern int  PPIDD_Uses_ma();
 extern int  PPIDD_MA_init(int dtype, int64_t *stack, int64_t *heap);
 extern double PPIDD_Wtime();
 extern void PPIDD_Error(char *message, int code);
 extern void PPIDD_Helper_server(int flag, int num_per_server);
 extern void PPIDD_Nxtval(int numproc, int64_t *val);
 extern void PPIDD_Size_all(int64_t *np);
 extern void PPIDD_Size(int64_t *np);
 extern void PPIDD_Rank(int64_t *me);
 extern void PPIDD_Init_fence();
 extern void PPIDD_Fence();
 extern void PPIDD_Send(void *buf,int64_t *count,int dtype,int64_t *dest,int64_t *sync);
 extern void PPIDD_Recv(void *buf,int64_t *count,int dtype,int64_t *source,int64_t *lenreal,int64_t *sourcereal,int64_t *sync);
 extern void PPIDD_Wait(int64_t *nodesel);
 extern int  PPIDD_Iprobe(int64_t *tag,int64_t *source);
 extern void PPIDD_BCast(void *buffer,int64_t *count,int dtype,int64_t *root);
 extern void PPIDD_Barrier();
 extern void PPIDD_Gsum(int dtype,void *buffer,int64_t *len, char *op);
 extern int  PPIDD_Create_irreg(char *name,int64_t *lenin,int64_t *nchunk,int dtype,int64_t *storetype,int64_t *handle);
 extern int  PPIDD_Create(char *name,int64_t *lentot, int dtype, int64_t *storetype, int64_t *handle);
 extern int  PPIDD_Destroy(int64_t *handle);
 extern int  PPIDD_Distrib(int64_t *handle,int64_t *rank,int64_t *ilo,int64_t *ihi);
 extern int  PPIDD_Location(int64_t *handle,int64_t *ilo,int64_t *ihi,int64_t *map,int64_t *proclist,int64_t *npk);
 extern int  PPIDD_Get(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff);
 extern int  PPIDD_Put(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff);
 extern int  PPIDD_Acc(int64_t *handle,int64_t *ilo,int64_t *ihi,void *buff,void *fac);
 extern void PPIDD_Read_inc(int64_t *handle,int64_t *inum,int64_t *incr,int64_t *returnval);
 extern void PPIDD_Zero_patch(int64_t *handle,int64_t *ilo,int64_t *ihi);
 extern int  PPIDD_Zero(int64_t *handle);
 extern void PPIDD_Duplicate(int64_t *handlei, int64_t *handlej, char *name);
 extern void PPIDD_Inquire_name(int64_t *handle, char *name);
 extern void PPIDD_Inquire_stype(int64_t *handle, int64_t *storetype);
 extern void PPIDD_Inquire_mem(int64_t *mem_used);
 extern int  PPIDD_Create_mutexes(int64_t *storetype,int64_t *number);
 extern void PPIDD_Lock_mutex(int64_t *inum);
 extern void PPIDD_Unlock_mutex(int64_t *inum);
 extern int  PPIDD_Destroy_mutexes();

 extern int  PPIDD_Eaf_open(char *name, int64_t *type, int64_t *handle);
 extern int  PPIDD_Eaf_write(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length);
 extern int  PPIDD_Eaf_awrite(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id);
 extern int  PPIDD_Eaf_read(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length);
 extern int  PPIDD_Eaf_aread(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id);
 extern int  PPIDD_Eaf_wait(int64_t *handle,int64_t *request_id);
 extern int  PPIDD_Eaf_waitall(int64_t *list, int64_t *num);
 extern int  PPIDD_Eaf_probe(int64_t *request_id,int64_t *status);
 extern int  PPIDD_Eaf_close(int64_t *handle);
 extern int  PPIDD_Eaf_delete(char *name);
 extern int  PPIDD_Eaf_length(int64_t *handle,double *fsize);
 extern int  PPIDD_Eaf_truncate(int64_t *handle,double *offset);
 extern void PPIDD_Eaf_errmsg(int *code,char *message);

 extern int  PPIDD_Sf_create(char *name, double *size_hard_limit, double *size_soft_limit, double *req_size, int64_t *handle);
 extern int  PPIDD_Sf_write(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id);
 extern int  PPIDD_Sf_read(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id);
 extern int  PPIDD_Sf_wait(int64_t *request_id);
 extern int  PPIDD_Sf_waitall(int64_t *list, int64_t *num);
 extern int  PPIDD_Sf_destroy(int64_t *handle);
 extern void PPIDD_Sf_errmsg(int *code,char *message);