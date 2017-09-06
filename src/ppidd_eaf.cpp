#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

#ifdef MPI2
#include <mpi.h>
#define   MPI_EAF_RW -1
#define   MPI_EAF_W  -2
#define   MPI_EAF_R  -3
extern MPI_Comm MPIGA_WORK_COMM;
#elif defined(GA_MPI)
#include <ga.h>
extern "C" {
#include <eaf.h>
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
   Creates an EAF file using name and path specified in name as a template.
   Return the EAF file descriptor in handle.
   It is a non-collective operation.
\* ************************************************************************ */
   int PPIDD_Eaf_open(char *name,int64_t *type, int64_t *handle) {
#if defined(MPI2) || defined(GA_MPI)
#ifdef MPI2
      MPI_File mpi_fh;
      MPI_Fint mpifhandle;
      MPI_Comm mpicomm=MPIGA_WORK_COMM;
      int amode=0;
#endif
#ifdef GA_MPI
      int gahandle;
#endif
      int modetype=(int)*type;

      if(MPI_Debug)printf("In PPIDD_Eaf_open: begin.\n");
#ifdef MPI2
      switch(modetype){
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
      *handle=(int64_t)mpifhandle;
#endif
#ifdef GA_MPI
      int ierr=EAF_Open(name, modetype, &gahandle);
      *handle=(int64_t)gahandle;
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_open: end. handle=%d,ierr=%d\n",(int)*handle,ierr);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_open is not available in serial case.\n");
      exit(1);
#endif
   }

/* ******************************************************************************************** *\
   Synchronously write to the file specified by the file handle.
   Writes number of bytes to the file identified by handle at location offset.
\* ******************************************************************************************** */
   int PPIDD_Eaf_write(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,*byte_offset,(long)*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)(*byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : before MPI_File_write_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_write_at(mpi_fh,offset,buff,count,datatype,&status);
      if(MPI_Debug)printf("In PPIDD_Eaf_write : end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;

      int ierr=EAF_Write(gahandle,offset,buff,bytes);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_write is not available in serial case.\n");
      exit(1);
#endif
   }

/* ******************************************************************************************** *\
   Asynchronously write to the file specified by the file handle, and return a handle to the asynchronous operation.
   Writes number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when eaf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   int PPIDD_Eaf_awrite(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request;
#else
      MPIO_Request request;
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,*byte_offset,(long)*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)(*byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : before MPI_File_iwrite_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_iwrite_at(mpi_fh,offset,buff,count,datatype,&request);
      *request_id=(int64_t)request;
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,ierr,(int)*request_id,(long)request);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int request;

      int ierr=EAF_Awrite(gahandle,offset,buff,bytes,&request);
      *request_id=(int64_t)request;
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_awrite is not available in serial case.\n");
      exit(1);
#endif
   }


/* ******************************************************************************************** *\
   Synchronously read from the file specified by the file handle.
   Reads number of bytes to the file identified by handle at location offset.
\* ******************************************************************************************** */
   int PPIDD_Eaf_read(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,*byte_offset,(long)*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)(*byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : before MPI_File_read_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_read_at(mpi_fh,offset,buff,count,datatype,&status);
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;

      int ierr=EAF_Read(gahandle,offset,buff,bytes);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_read is not available in serial case.\n");
      exit(1);
#endif
   }


/* ******************************************************************************************** *\
   Asynchronously read from the file specified by the file handle, and return a handle to the asynchronous operation.
   Reads number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when eaf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   int PPIDD_Eaf_aread(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      MPI_Datatype datatype;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request;
#else
      MPIO_Request request;
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,*byte_offset,(long)*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)(*byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : before MPI_File_iread_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      int ierr=MPI_File_iread_at(mpi_fh,offset,buff,count,datatype,&request);
      *request_id=(int64_t)request;
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,ierr,(int)*request_id,(long)request);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int request;

      int ierr=EAF_Aread(gahandle,offset,buff,bytes,&request);
      *request_id=(int64_t)request;
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_aread is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************************************************ *\
   Wait for the completion of the asynchronous request associated with request_id.
   Request_id is invalidated after the calling.
   integer request_id   --[in]  Handle of asynchronous request.
   integer ierr         --[out] Error code. 0 if it is able to wait for completion,
                          else returns error code.
\* ************************************************************************************ */
   int PPIDD_Eaf_wait(int64_t *handle,int64_t *request_id) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request=(MPI_Request)(*request_id);
#else
      MPIO_Request request=(MPIO_Request)(*request_id);
#endif
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_wait  : begin. handle=%d,request_id=%d,request=%ld\n",(int)mpifhandle,(int)*request_id,(long)request);

#ifdef MPIO_USES_MPI_REQUEST
      int ierr=MPI_Wait( &request, &status );
#else
      int ierr=MPIO_Wait(&request, &status);
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_wait  : end. ierr=%d\n",ierr);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      int request=(int)*request_id;

      int ierr=EAF_Wait(gahandle,request);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_wait is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************** *\
   Blocks the calling process until all of the num I/O operations associated with ids
   specified in list complete. Finally invalidates (modifies) ids on the list.
\* ********************************************************************************** */
   int PPIDD_Eaf_waitall(int64_t *list, int64_t *num) {
#ifdef MPI2
      int count=(int)(*num);
      MPI_Status *array_of_statuses;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request *array_of_requests;
      array_of_requests=(MPI_Request*)malloc(count*sizeof(MPI_Request));
#else
      MPIO_Request *array_of_requests;
      array_of_requests=(MPIO_Request*)malloc(count*sizeof(MPIO_Request));
#endif

      array_of_statuses=(MPI_Status*)malloc(count*sizeof(MPI_Status));
      for(int i=0;i<count;i++) array_of_requests[i]=(MPI_Request)list[i];

#ifdef MPIO_USES_MPI_REQUEST
      int ierr=MPI_Waitall(count,array_of_requests,array_of_statuses);
#else
      int ierr=0;
      for(int i=0;i<count;i++) ierr+=MPIO_Wait(&array_of_requests[i], &array_of_statuses[i]);
#endif
      return ierr;
#elif defined(GA_MPI)
      return 0;
#else
      printf(" ERROR: PPIDD_Eaf_waitall is not available in serial case.\n");
      exit(1);
#endif
   }

/* ************************************************************************************ *\
   Determine if an asynchronous request has completed or is pending.
   integer request_id   --[in]  Handle of asynchronous request.
   integer status       --[out] Pending or completed status argument.
                          status returns 0 if the asynchronous operation is complete, or 1 otherwise.
                          If the asynchronous request is complete, id is invalidated.
   integer ierr         --[out] Error code. 0 if probe succeeded, else returns error code.
\* ************************************************************************************ */
   int PPIDD_Eaf_probe(int64_t *request_id,int64_t *status) {
#ifdef MPI2
      int flag;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request=(MPI_Request)(*request_id);
#else
      MPIO_Request request=(MPIO_Request)(*request_id);
#endif
      MPI_Status mpistatus;
      if(MPI_Debug)printf("In PPIDD_Eaf_probe  : begin. request_id=%d,request=%ld\n",(int)*request_id,(long)request);

#ifdef MPIO_USES_MPI_REQUEST
      int ierr=MPI_Test(&request, &flag, &mpistatus);
#else
      int ierr=MPIO_Test(&request, &flag, &mpistatus);
#endif
      if(flag) *status=(int64_t)0;
      else *status=(int64_t)1;

      if(MPI_Debug)printf("In PPIDD_Eaf_probe  : end. ierr=%d\n",ierr);
      return ierr;
#elif defined(GA_MPI)
      int garequest=(int)*request_id;
      int gastatus;

      int ierr=EAF_Probe(garequest, &gastatus);
      *status=(int64_t)gastatus;
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_probe is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************************************************ *\
   Close an EAF file.
   integer handle  --[in]  File Handle.
   integer ierr    --[out] Error code. 0 if the file was closed, else returns error code.
\* ************************************************************************************ */
   int PPIDD_Eaf_close(int64_t *handle) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      if(MPI_Debug)printf("In PPIDD_Eaf_close: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_close( &mpi_fh );
      if(MPI_Debug)printf("In PPIDD_Eaf_close: end. handle=%d,ierr=%d\n",(int)mpifhandle,ierr);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;

      int ierr=EAF_Close(gahandle);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_close is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************************* *\
   Delete an EAF file
   character*(*) name   -- [in]  File name.
   integer ierr         -- [out] Error code. 0 if the file was deleted, else returns error code.
\* ********************************************************************************************* */
   int PPIDD_Eaf_delete(char *name) {
#if defined(MPI2) || defined(GA_MPI)
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: begin. name=%s\n",name);
#ifdef MPI2
      int ierr=MPI_File_delete(name,MPI_INFO_NULL);
#elif defined(GA_MPI)
      int ierr=EAF_Delete(name);
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: end. ierr=%d\n",ierr);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_delete is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************************************************ *\
   Determine the length (in bytes) of an EAF file.
   integer handle    --[in]  File Handle.
   double fsize      --[out] File length in bytes.
\* ************************************************************************************ */
   int PPIDD_Eaf_length(int64_t *handle,double *fsize) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset size;
      if(MPI_Debug)printf("In PPIDD_Eaf_length: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_get_size(mpi_fh, &size);
      *fsize=(double)size;
      if(MPI_Debug)printf("In PPIDD_Eaf_length: end. handle=%d,fsize=%f\n",(int)mpifhandle,*fsize);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t length;

      int ierr=EAF_Length(gahandle,&length);
      *fsize=(double)length;
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_length is not available in serial case.\n");
      exit(1);
#endif
   }

/* *************************************************************************************** *\
   Truncate an EAF file at specified offset (in bytes).
   integer handle --[in]  File Handle.
   double offset  --[in]  Offset in bytes.
   integer ierr   --[out] Error code. 0 if the file was truncated, else returns error code.
\* *************************************************************************************** */
   int PPIDD_Eaf_truncate(int64_t *handle,double *offset) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset size=(MPI_Offset)(*offset);
      if(MPI_Debug)printf("In PPIDD_Eaf_truncate: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      int ierr=MPI_File_set_size(mpi_fh,size);
      if(MPI_Debug)printf("In PPIDD_Eaf_truncate: end. handle=%d,size=%ld\n",(int)mpifhandle,(long)size);
      return ierr;
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t length=(eaf_off_t)*offset;

      int ierr=EAF_Truncate(gahandle,length);
      return ierr;
#else
      printf(" ERROR: PPIDD_Eaf_truncate is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************* *\
   Returns a string interpretation of the error code, or an empty string
   (Fortran all blanks, C null first character) if the error code is not recognized.
        code             -- [in]  error code returned by a previous call to EAF
        message          -- [out] character string where the corresponding message
\* ********************************************************************************* */
   void PPIDD_Eaf_errmsg(int *code,char *message) {
#if defined(MPI2) || defined(GA_MPI)
#ifdef MPI2
      int eclass, len;
      char estring[MPI_MAX_ERROR_STRING],estring2[MPI_MAX_ERROR_STRING];
#endif
      int lxi=strlen(message);

      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: begin. code=%d\n",*code);
#ifdef MPI2
      MPI_Error_class(*code, &eclass);
      MPI_Error_string(*code, estring, &len);
      sprintf(estring2," Error %d: %s", eclass, estring);
      strcpy(message,estring2);
#elif defined(GA_MPI)
      EAF_Errmsg(*code, message);
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: middle. message=%s\n",message);
      for(int i=strlen(message);i<lxi;i++) message[i]=' ';
      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: end. code=%d\n",*code);
#else
      printf(" ERROR: PPIDD_Eaf_errmsg is not available in serial case.\n");
      exit(1);
#endif
   }
