#ifndef __PPIDD_SF_H__
#define __PPIDD_SF_H__
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
 extern void PPIDD_Sf_create(char *fname, double *size_hard_limit, double *size_soft_limit, double *req_size, int64_t *handle, int64_t *ierr);
 extern void PPIDD_Sf_write(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id,int64_t *ierr);
 extern void PPIDD_Sf_read(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id,int64_t *ierr);
 extern void PPIDD_Sf_wait(int64_t *request_id,int64_t *ierr);
 extern void PPIDD_Sf_waitall(int64_t *list, int64_t *num,int64_t *ierr);
 extern void PPIDD_Sf_destroy(int64_t *handle,int64_t *ierr);
 extern void PPIDD_Sf_errmsg(int64_t *ierr,char *message);
#ifdef __cplusplus
}
#endif
#endif
