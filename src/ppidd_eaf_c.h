#ifndef __PPIDD_EAF_C_H__
#define __PPIDD_EAF_C_H__
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
 extern void PPIDD_Eaf_open(char *fname, int64_t *type, int64_t *handle, int64_t *ierr);
 extern void PPIDD_Eaf_write(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *ierr);
 extern void PPIDD_Eaf_awrite(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id,int64_t *ierr);
 extern void PPIDD_Eaf_read(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *ierr);
 extern void PPIDD_Eaf_aread(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id,int64_t *ierr);
 extern void PPIDD_Eaf_wait(int64_t *handle,int64_t *request_id,int64_t *ierr);
 extern void PPIDD_Eaf_waitall(int64_t *list, int64_t *num,int64_t *ierr);
 extern void PPIDD_Eaf_probe(int64_t *request_id,int64_t *status,int64_t *ierr);
 extern void PPIDD_Eaf_close(int64_t *handle,int64_t *ierr);
 extern void PPIDD_Eaf_delete(char *fname, int64_t *ierr);
 extern void PPIDD_Eaf_length(int64_t *handle,double *fsize,int64_t *ierr);
 extern void PPIDD_Eaf_truncate(int64_t *handle,double *offset,int64_t *ierr);
 extern void PPIDD_Eaf_errmsg(int64_t *code,char *message);
#ifdef __cplusplus
}
#endif
#endif
