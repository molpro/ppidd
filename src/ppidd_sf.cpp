#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

#ifdef MPI2
#include <mpi.h>
extern MPI_Comm MPIGA_WORK_COMM;
#elif defined(GA_MPI)
#include <ga.h>
extern "C" {
#include <sf.h>
}
#endif

#include "ppidd.h"

#if defined(GA_MPI) || defined(MPI2)
static int MPI_Debug=0;
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ************************************************************************ *\
   Creates shared file using name and path specified in fname as a template.
   req_size specifies size of a typical request (-1. means "don't know").
   It is a collective operation.
\* ************************************************************************ */
   void PPIDD_Sf_create(char *fname ,double *size_hard_limit, double *size_soft_limit, double *req_size, int64_t *handle, int64_t *ierr) {
#if defined(MPI2) || defined(GA_MPI)
#ifdef MPI2
      MPI_Comm mpicomm=MPIGA_WORK_COMM;
      MPI_File mpi_fh;
      MPI_Fint mpifhandle;
      int mpierr;
#elif defined(GA_MPI)
      char *errmsg;
#endif
      int i;
      char *name2;
      int lxi;

      if(MPI_Debug)printf("In PPIDD_Sf_create: begin.\n");
      lxi=strlen(fname);
      strncpy((name2=(char *)malloc(lxi+1)),fname,lxi);
      name2[lxi]=(char)0;
      for(i=lxi-1; (i>=0 && name2[i]==' '); i--) name2[i]=(char)0;

#ifdef MPI2
      mpierr=MPI_File_open(mpicomm,name2,MPI_MODE_RDWR|MPI_MODE_CREATE|MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_UNIQUE_OPEN,MPI_INFO_NULL,&mpi_fh);
      mpifhandle = MPI_File_c2f( mpi_fh );
      *handle=(int64_t)mpifhandle;
      *ierr=(int64_t)mpierr;
#elif defined(GA_MPI)
      if(MPI_Debug) {
         printf("PPIDD_Sf_create: sizeof(double) =%d, sizeof(SFsize_t)=%d\n",(int)sizeof(double),(int)sizeof(SFsize_t));
         printf("PPIDD_Sf_create: sizeof(int64_t)=%d, sizeof(Integer) =%d\n",(int)sizeof(int64_t),(int)sizeof(Integer));
      }
      if ( sizeof(double) != sizeof(SFsize_t) ) {
         printf("PPIDD_Sf_create: sizeof(double) =%d, sizeof(SFsize_t)=%d\n",(int)sizeof(double),(int)sizeof(SFsize_t));
         errmsg=strdup(" PPIDD_Sf_create: Data types do not match between [double] and [SFsize_t]");
         GA_Error(errmsg,0);
         free(errmsg);
      }
      if ( sizeof(int64_t) != sizeof(Integer) ) {
         printf("PPIDD_Sf_create: sizeof(int64_t)=%d, sizeof(Integer) =%d\n",(int)sizeof(int64_t),(int)sizeof(Integer));
         errmsg=strdup(" PPIDD_Sf_create: Data types do not match between [int64_t] and [Integer]");
         GA_Error(errmsg,0);
         free(errmsg);
      }

      *ierr=(int64_t)SF_Create(name2, *size_hard_limit, *size_soft_limit, *req_size, &i);
      *handle=(int64_t)i;
#endif

      free(name2);
      if(MPI_Debug)printf("In PPIDD_Sf_create: end. handle=%d,ierr=%d\n",(int)*handle,(int)*ierr);
#else
      printf(" ERROR: PPIDD_Sf_create is not available in serial case.\n");
      exit(1);
#endif
   }


/* ******************************************************************************************** *\
   Asynchronous write operation.
   Writes number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when sf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   void PPIDD_Sf_write(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      int mpierr;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request;
#else
      MPIO_Request request;
#endif
      if(MPI_Debug)printf("In PPIDD_Sf_write : begin. handle=%d,byte_offset=%f,byte_length=%f\n",(int)mpifhandle,*byte_offset,*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)((*byte_length)/sizeof(double));
      if(MPI_Debug)printf("In PPIDD_Sf_write : before MPI_File_iwrite_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      mpierr=MPI_File_iwrite_at(mpi_fh,offset,buff,count,MPI_DOUBLE,&request);
      *request_id=(int64_t)request;
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Sf_write : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,(int)*ierr,(int)*request_id,(long)request);
#elif defined(GA_MPI)
      char *buffer=(char *)buff;

      int ihandle=(int)*handle;
      int irequest_id;
      *ierr=(int64_t)SF_Write(ihandle, *byte_offset, *byte_length, buffer, &irequest_id);
      *request_id=(int64_t)irequest_id;
#else
      printf(" ERROR: PPIDD_Sf_write is not available in serial case.\n");
      exit(1);
#endif
   }

/* ******************************************************************************************** *\
   Asynchronous read operation.
   Reads number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when sf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   void PPIDD_Sf_read(int64_t *handle,double *byte_offset,double *byte_length, double *buff,int64_t *request_id,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      int mpierr;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request;
#else
      MPIO_Request request;
#endif
      if(MPI_Debug) printf("In PPIDD_Sf_read  : begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)((*byte_length)/sizeof(double));
      if(MPI_Debug)printf("In PPIDD_Sf_read  : before MPI_File_iread_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      mpierr=MPI_File_iread_at(mpi_fh,offset,buff,count,MPI_DOUBLE,&request);
      *request_id=(int64_t)request;
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Sf_read  : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,(int)*ierr,(int)*request_id,(long)request);
#elif defined(GA_MPI)
      char *buffer=(char *)buff;

      int ihandle=(int)*handle;
      int irequest_id;
      *ierr=(int64_t)SF_Read(ihandle, *byte_offset, *byte_length, buffer, &irequest_id);
      *request_id=(int64_t)irequest_id;
#else
      printf(" ERROR: PPIDD_Sf_read is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************************************************ *\
   Blocks the calling process until I/O operation associated with request_id completes.
   Invalidates request_id.
\* ************************************************************************************ */
   void PPIDD_Sf_wait(int64_t *request_id,int64_t *ierr) {
#ifdef MPI2
      int mpierr;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request=(MPI_Request)(*request_id);
#else
      MPIO_Request request=(MPIO_Request)(*request_id);
#endif
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Sf_wait  : begin. request_id=%d,request=%ld\n",(int)*request_id,(long)request);
#ifdef MPIO_USES_MPI_REQUEST
      mpierr=MPI_Wait( &request, &status );
#else
      mpierr=MPIO_Wait(&request, &status);
#endif
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Sf_wait  : end. ierr=%d\n",(int)*ierr);
#elif defined(GA_MPI)
      int irequest_id=(int)*request_id;
      *ierr=(int64_t)SF_Wait(&irequest_id);
      *request_id=(int64_t)irequest_id;
#else
      printf(" ERROR: PPIDD_Sf_wait is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************** *\
   Blocks the calling process until all of the num I/O operations associated with ids
   specified in list complete. Invalidates (modifies) ids on the list.
\* ********************************************************************************** */
   void PPIDD_Sf_waitall(int64_t *list, int64_t *num,int64_t *ierr) {
#ifdef MPI2
      int i;
      int mpierr;
      int count=(int)(*num);
      MPI_Status *array_of_statuses;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request *array_of_requests;
      array_of_requests=(MPI_Request*)malloc(count*sizeof(MPI_Request));
#else
      int mpierrsub;
      MPIO_Request *array_of_requests;
      array_of_requests=(MPIO_Request*)malloc(count*sizeof(MPIO_Request));
#endif

      array_of_statuses=(MPI_Status*)malloc(count*sizeof(MPI_Status));
      for(i=0;i<count;i++) array_of_requests[i]=(MPI_Request)list[i];

#ifdef MPIO_USES_MPI_REQUEST
      mpierr=MPI_Waitall(count,array_of_requests,array_of_statuses);
#else
      for(i=0,mpierr=0;i<count;i++) {
        mpierrsub=MPIO_Wait(&array_of_requests[i], &array_of_statuses[i]);
        mpierr=mpierr+mpierrsub;
      }
#endif
      *ierr=(int64_t)mpierr;
#elif defined(GA_MPI)
      int inum=(int)*num;
      int *ilist,i;

      ilist=(int *)malloc(inum * sizeof(int));
      for (i=0; i<inum; ++i) {
          ilist[i] = (int) list[i];
      }

      *ierr=(int64_t)SF_Waitall(ilist, inum);

      for (i=0; i<inum; ++i) {
          list[i] = (int64_t) ilist[i];
      }
      free(ilist);
#else
      printf(" ERROR: PPIDD_Sf_waitall is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************ *\
   Destroys the shared file associated with handle.
   It is a collective operation.
\* ************************************************ */
   void PPIDD_Sf_destroy(int64_t *handle,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      int mpierr;
      if(MPI_Debug)printf("In PPIDD_Sf_destroy: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      mpierr=MPI_File_close( &mpi_fh );
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Sf_destroy: end. handle=%d,ierr=%d\n",(int)mpifhandle,(int)*ierr);
#elif defined(GA_MPI)
      int ihandle=(int)*handle;
      *ierr=(int64_t)SF_Destroy(ihandle);
#else
      printf(" ERROR: PPIDD_Sf_destroy is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************* *\
   Returns a string interpretation of the error code, or an empty string
   (Fortran all blanks, C null first character) if the error code is not recognized.
        ierr             -- error code returned by a previous call to SF [in]
        message          -- character string where the corresponding message
                            is written [out]
\* ********************************************************************************* */
   void PPIDD_Sf_errmsg(int64_t *ierr,char *message) {
#if defined(MPI2) || defined(GA_MPI)
#ifdef MPI2
      int eclass, len;
      char estring[MPI_MAX_ERROR_STRING],estring2[MPI_MAX_ERROR_STRING];
#endif
      int perrcode=(int)*ierr;
      int i;
      int lxi;

      lxi=strlen(message);

      if(MPI_Debug)printf("In PPIDD_Sf_errmsg: begin. perrcode=%d\n",perrcode);
#ifdef MPI2
      MPI_Error_class(perrcode, &eclass);
      MPI_Error_string(perrcode, estring, &len);
      sprintf(estring2," Error %d: %s", eclass, estring);
      strcpy(message,estring2);
#elif defined(GA_MPI)
      SF_Errmsg(perrcode, message);
#endif
      if(MPI_Debug)printf("In PPIDD_Sf_errmsg: middle. message=%s\n",message);
      for(i=strlen(message);i<lxi;i++) message[i]=' ';
      if(MPI_Debug)printf("In PPIDD_Sf_errmsg: end. perrcode=%d\n",perrcode);
#else
      printf(" ERROR: PPIDD_Sf_errmsg is not available in serial case.\n");
      exit(1);
#endif
   }