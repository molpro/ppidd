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
   Creates an EAF file using name and path specified in fname as a template.
   Return the EAF file descriptor in handle.
   It is a non-collective operation.
\* ************************************************************************ */
   void PPIDD_Eaf_open(char *fname,int64_t *type, int64_t *handle, int64_t *ierr) {
#if defined(MPI2) || defined(GA_MPI)
#ifdef MPI2
      MPI_File mpi_fh;
      MPI_Fint mpifhandle;
      MPI_Comm mpicomm=MPIGA_WORK_COMM;
      int amode=0;
      int mpierr;
#endif
#ifdef GA_MPI
      int gaerr;
      int gahandle;
#endif
      int modetype=(int)*type;
      int i;
      char *name2;
      int lxi;


      if(MPI_Debug)printf("In PPIDD_Eaf_open: begin.\n");
      lxi=strlen(fname);
      strncpy((name2=(char *)malloc(lxi+1)),fname,lxi);
      name2[lxi]=(char)0;
      for(i=lxi-1; (i>=0 && name2[i]==' '); i--) name2[i]=(char)0;

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

      mpierr=MPI_File_open(MPI_COMM_SELF,name2,amode,MPI_INFO_NULL,&mpi_fh);
      mpifhandle = MPI_File_c2f( mpi_fh );
      *handle=(int64_t)mpifhandle;
      *ierr=(int64_t)mpierr;
#endif
#ifdef GA_MPI
      gaerr=EAF_Open(name2, modetype, &gahandle);
      *handle=(int64_t)gahandle;
      *ierr=(int64_t)gaerr;
#endif
      free(name2);
      if(MPI_Debug)printf("In PPIDD_Eaf_open: end. handle=%d,ierr=%d\n",(int)*handle,(int)*ierr);
#else
      printf(" ERROR: PPIDD_Eaf_open is not available in serial case.\n");
      exit(1);
#endif
   }

/* ******************************************************************************************** *\
   Synchronously write to the file specified by the file handle.
   Writes number of bytes to the file identified by handle at location offset.
\* ******************************************************************************************** */
   void PPIDD_Eaf_write(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      int mpierr;
      MPI_Datatype datatype;
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,*byte_offset,(long)*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)(*byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : before MPI_File_write_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      mpierr=MPI_File_write_at(mpi_fh,offset,buff,count,datatype,&status);
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_write : end. handle=%d,ierr=%d\n",(int)mpifhandle,(int)*ierr);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int gaerr;

      gaerr=EAF_Write(gahandle,offset,buff,bytes);
      *ierr=(int64_t)gaerr;
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
   void PPIDD_Eaf_awrite(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      int mpierr;
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
      mpierr=MPI_File_iwrite_at(mpi_fh,offset,buff,count,datatype,&request);
      *request_id=(int64_t)request;
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_awrite : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,(int)*ierr,(int)*request_id,(long)request);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int request;
      int gaerr;

      gaerr=EAF_Awrite(gahandle,offset,buff,bytes,&request);
      *request_id=(int64_t)request;
      *ierr=(int64_t)gaerr;
#else
      printf(" ERROR: PPIDD_Eaf_awrite is not available in serial case.\n");
      exit(1);
#endif
   }


/* ******************************************************************************************** *\
   Synchronously read from the file specified by the file handle.
   Reads number of bytes to the file identified by handle at location offset.
\* ******************************************************************************************** */
   void PPIDD_Eaf_read(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      int mpierr;
      MPI_Datatype datatype;
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : begin. handle=%d,byte_offset=%f,byte_length=%ld\n",(int)mpifhandle,*byte_offset,(long)*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)(*byte_length/8);
      datatype=MPI_DOUBLE;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : before MPI_File_read_at. handle=%d,offset=%ld,count=%d\n",(int)mpifhandle,(long)offset,count);
      mpierr=MPI_File_read_at(mpi_fh,offset,buff,count,datatype,&status);
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_read  : end. handle=%d,ierr=%d\n",(int)mpifhandle,(int)*ierr);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int gaerr;

      gaerr=EAF_Read(gahandle,offset,buff,bytes);
      *ierr=(int64_t)gaerr;
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
   void PPIDD_Eaf_aread(int64_t *handle,double *byte_offset,void *buff,int64_t *byte_length,int64_t *request_id,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset offset;
      int count;
      int mpierr;
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
      mpierr=MPI_File_iread_at(mpi_fh,offset,buff,count,datatype,&request);
      *request_id=(int64_t)request;
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_aread  : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",(int)mpifhandle,(int)*ierr,(int)*request_id,(long)request);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t offset=(eaf_off_t)*byte_offset;
      size_t bytes=(size_t)*byte_length;
      int request;
      int gaerr;

      gaerr=EAF_Aread(gahandle,offset,buff,bytes,&request);
      *request_id=(int64_t)request;
      *ierr=(int64_t)gaerr;
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
   void PPIDD_Eaf_wait(int64_t *handle,int64_t *request_id,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      int mpierr;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request=(MPI_Request)(*request_id);
#else
      MPIO_Request request=(MPIO_Request)(*request_id);
#endif
      MPI_Status status;
      if(MPI_Debug)printf("In PPIDD_Eaf_wait  : begin. handle=%d,request_id=%d,request=%ld\n",(int)mpifhandle,(int)*request_id,(long)request);

#ifdef MPIO_USES_MPI_REQUEST
      mpierr=MPI_Wait( &request, &status );
#else
      mpierr=MPIO_Wait(&request, &status);
#endif
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_wait  : end. ierr=%d\n",(int)*ierr);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      int request=(int)*request_id;
      int gaerr;

      gaerr=EAF_Wait(gahandle,request);
      *ierr=(int64_t)gaerr;
#else
      printf(" ERROR: PPIDD_Eaf_wait is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************** *\
   Blocks the calling process until all of the num I/O operations associated with ids
   specified in list complete. Finally invalidates (modifies) ids on the list.
\* ********************************************************************************** */
   void PPIDD_Eaf_waitall(int64_t *list, int64_t *num,int64_t *ierr) {
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
   void PPIDD_Eaf_probe(int64_t *request_id,int64_t *status,int64_t *ierr) {
#ifdef MPI2
      int flag;
      int mpierr;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request=(MPI_Request)(*request_id);
#else
      MPIO_Request request=(MPIO_Request)(*request_id);
#endif
      MPI_Status mpistatus;
      if(MPI_Debug)printf("In PPIDD_Eaf_probe  : begin. request_id=%d,request=%ld\n",(int)*request_id,(long)request);

#ifdef MPIO_USES_MPI_REQUEST
      mpierr=MPI_Test(&request, &flag, &mpistatus);
#else
      mpierr=MPIO_Test(&request, &flag, &mpistatus);
#endif
      if(flag) *status=(int64_t)0;
      else *status=(int64_t)1;

      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_probe  : end. ierr=%d\n",(int)*ierr);
#elif defined(GA_MPI)
      int garequest=(int)*request_id;
      int gastatus;
      int gaerr;

      gaerr=EAF_Probe(garequest, &gastatus);
      *status=(int64_t)gastatus;
      *ierr=(int64_t)gaerr;
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
   void PPIDD_Eaf_close(int64_t *handle,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      int mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_close: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      mpierr=MPI_File_close( &mpi_fh );
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_close: end. handle=%d,ierr=%d\n",(int)mpifhandle,(int)*ierr);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      int gaerr;

      gaerr=EAF_Close(gahandle);
      *ierr=(int64_t)gaerr;
#else
      printf(" ERROR: PPIDD_Eaf_close is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************************* *\
   Delete an EAF file
   character*(*) fname   -- [in]  File name.
   integer ierr          -- [out] Error code. 0 if the file was deleted, else returns error code.
\* ********************************************************************************************* */
   void PPIDD_Eaf_delete(char *fname,int64_t *ierr) {
#if defined(MPI2) || defined(GA_MPI)
      int i,pierr=0;
      char *name2;
      int lxi;

      if(MPI_Debug)printf("In PPIDD_Eaf_delete: begin.\n");
      lxi=strlen(fname);
      strncpy((name2=(char *)malloc(lxi+1)),fname,lxi);
      name2[lxi]=(char)0;
      for(i=lxi-1; (i>=0 && name2[i]==' '); i--) name2[i]=(char)0;

#ifdef MPI2
      pierr=MPI_File_delete(name2,MPI_INFO_NULL);
#elif defined(GA_MPI)
      pierr=EAF_Delete(name2);
#endif

      if(MPI_Debug)printf("In PPIDD_Eaf_delete: mid. fname=%s,pierr=%d\n",name2,pierr);
      free(name2);
      *ierr=(int64_t)pierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_delete: end. ierr=%d\n",(int)*ierr);
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
   void PPIDD_Eaf_length(int64_t *handle,double *fsize,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset size;
      int mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_length: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      mpierr=MPI_File_get_size(mpi_fh, &size);
      *fsize=(double)size;
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_length: end. handle=%d,fsize=%f\n",(int)mpifhandle,*fsize);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t length;
      int gaerr;

      gaerr=EAF_Length(gahandle,&length);
      *fsize=(double)length;
      *ierr=(int64_t)gaerr;
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
   void PPIDD_Eaf_truncate(int64_t *handle,double *offset,int64_t *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      MPI_Offset size=(MPI_Offset)(*offset);
      int mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_truncate: begin. handle=%d\n",(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      mpierr=MPI_File_set_size(mpi_fh,size);
      *ierr=(int64_t)mpierr;
      if(MPI_Debug)printf("In PPIDD_Eaf_truncate: end. handle=%d,size=%ld\n",(int)mpifhandle,(long)size);
#elif defined(GA_MPI)
      int gahandle=(int)*handle;
      eaf_off_t length=(eaf_off_t)*offset;
      int gaerr;

      gaerr=EAF_Truncate(gahandle,length);
      *ierr=(int64_t)gaerr;
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
   void PPIDD_Eaf_errmsg(int64_t *code,char *message) {
#if defined(MPI2) || defined(GA_MPI)
#ifdef MPI2
      int eclass, len;
      char estring[MPI_MAX_ERROR_STRING],estring2[MPI_MAX_ERROR_STRING];
#endif
      int perrcode=(int)*code;
      int i;
      int lxi;

      lxi=strlen(message);

      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: begin. perrcode=%d\n",perrcode);
#ifdef MPI2
      MPI_Error_class(perrcode, &eclass);
      MPI_Error_string(perrcode, estring, &len);
      sprintf(estring2," Error %d: %s", eclass, estring);
      strcpy(message,estring2);
#elif defined(GA_MPI)
      EAF_Errmsg(perrcode, message);
#endif
      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: middle. message=%s\n",message);
      for(i=strlen(message);i<lxi;i++) message[i]=' ';
      if(MPI_Debug)printf("In PPIDD_Eaf_errmsg: end. perrcode=%d\n",perrcode);
#else
      printf(" ERROR: PPIDD_Eaf_errmsg is not available in serial case.\n");
      exit(1);
#endif
   }
