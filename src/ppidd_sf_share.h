
/* ====================================================================== *\
 *                    PPIDD Shared Files Library                          *
 *                    ==========================                          *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ------------------------------------------------------------------------------------------ *
 * The Shared Files (SF) library implements logically-shared temporary files for parallel     *
 * SPMD (single-program-multiple-data) applications. The main features are listed as follows: *
 * -- Shared files are non-persistent (temporary)                                             *
 * -- Shared files resemble one-dimensional arrays in main memory                             *
 * -- Each process can independently read/write to any location in the file                   *
 * -- All routines return error code: "0" means success.                                      *
 * -- sf_create and sf_destroy are collective                                                 *
 * -- file, request sizes, and offset (all in bytes) are DOUBLE PRECISION arguments,          *
 *    all the other arguments are INTEGERS                                                    *
 * -- read/writes are asynchronous                                                            *
 *--------------------------------------------------------------------------------------------*
 * Common Fortran and C interfaces.                                       *
 * For Fortran interface, the subroutines in this file named PPIDD_Sf_XXX *
 * are converted to the proper FORTRAN external by the FORT_Extern macro  *
 * and the definitions in the ppidd_sf_fortran.h header file.             *
 * For C interface, the subroutines in this file named PPIDD_Sf_XXX can be*
 * only called by C program directly. Any calling by Fortran progam       *
 * should refer to the routines in the ppidd_sf_fortran.h file.           *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       15/07/2008                                                 *
\* ====================================================================== */



/* ************************************************************************ *\
   Creates shared file using name and path specified in fname as a template.
   req_size specifies size of a typical request (-1. means "don't know").
   It is a collective operation.
\* ************************************************************************ */
   void PPIDD_Sf_create(char *fname
#if defined(FORTCL_NEXT)
	       ,fortintc lx
#endif
	       ,double *size_hard_limit, double *size_soft_limit, double *req_size, fortint *handle, fortint *ierr
#if defined(FORTCL_END)
	       ,fortintc lx
#endif
)  {
#if defined(MPI2) || defined(GA_TOOLS)
#ifdef MPI2
      MPI_Comm mpicomm=MPIGA_WORK_COMM;
      MPI_File mpi_fh;
      MPI_Fint mpifhandle;
      int mpierr;
#endif
#ifdef GA_TOOLS
      char *errmsg;
#endif
      int i;
      char *name2;
      int lxi;

      if(MPI_Debug)printf("%5d: In PPIDD_Sf_create: begin.\n",ppidd_sf_rank());
#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(fname);
#endif
      strncpy((name2=(char *)malloc(lxi+1)),fname,lxi);
      name2[lxi]=(char)0;
      for(i=lxi-1; (i>=0 && name2[i]==' '); i--) name2[i]=(char)0;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_create: midlle.sizeof(fortint)=%d,sizeof(fortintc)=%d,lxi=%d,filename=%s\n",ppidd_sf_rank(),(int)sizeof(fortint),(int)sizeof(fortintc),lxi,name2);

#ifdef MPI2
      mpierr=MPI_File_open(mpicomm,name2,MPI_MODE_RDWR|MPI_MODE_CREATE|MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_UNIQUE_OPEN,MPI_INFO_NULL,&mpi_fh);
      mpifhandle = MPI_File_c2f( mpi_fh );
      *handle=(fortint)mpifhandle;
      *ierr=(fortint)mpierr;
#endif
#ifdef GA_TOOLS
      if(MPI_Debug) {
         printf("%5d: PPIDD_Sf_create: sizeof(double) =%d, sizeof(SFsize_t)=%d\n",ppidd_sf_rank(),(int)sizeof(double),(int)sizeof(SFsize_t));
         printf("%5d: PPIDD_Sf_create: sizeof(fortint)=%d, sizeof(Integer) =%d\n",ppidd_sf_rank(),(int)sizeof(fortint),(int)sizeof(Integer));
      }
      if ( sizeof(double) != sizeof(SFsize_t) ) {
         printf("%5d: PPIDD_Sf_create: sizeof(double) =%d, sizeof(SFsize_t)=%d\n",ppidd_sf_rank(),(int)sizeof(double),(int)sizeof(SFsize_t));
         errmsg=strdup(" PPIDD_Sf_create: Data types do not match between [double] and [SFsize_t]");
         GA_Error(errmsg,0);
         free(errmsg);
      }
      if ( sizeof(fortint) != sizeof(Integer) ) {
         printf("%5d: PPIDD_Sf_create: sizeof(fortint)=%d, sizeof(Integer) =%d\n",ppidd_sf_rank(),(int)sizeof(fortint),(int)sizeof(Integer));
         errmsg=strdup(" PPIDD_Sf_create: Data types do not match between [fortint] and [Integer]");
         GA_Error(errmsg,0);
         free(errmsg);
      }

      *ierr=(fortint)SF_Create(name2, *size_hard_limit, *size_soft_limit, *req_size, &i);
      *handle=(fortint)i;
#endif

      free(name2);
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_create: end. handle=%d,ierr=%d\n",ppidd_sf_rank(),(int)*handle,(int)*ierr);
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Sf_create is not available in serial case.\n");
      exit(1);
#endif
   }


/* ******************************************************************************************** *\
   Asynchronous write operation.
   Writes number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when sf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   void PPIDD_Sf_write(fortint *handle,double *byte_offset,double *byte_length, double *buff,fortint *request_id,fortint *ierr) {
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
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_write : begin. handle=%d,byte_offset=%f,byte_length=%f\n",ppidd_sf_rank(),(int)mpifhandle,*byte_offset,*byte_length);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)((*byte_length)/sizeof(double));
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_write : before MPI_File_iwrite_at. handle=%d,offset=%ld,count=%d\n",ppidd_sf_rank(),(int)mpifhandle,(long)offset,count);
      mpierr=MPI_File_iwrite_at(mpi_fh,offset,buff,count,MPI_DOUBLE,&request);
      *request_id=(fortint)request;
      *ierr=(fortint)mpierr;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_write : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",ppidd_sf_rank(),(int)mpifhandle,(int)*ierr,(int)*request_id,(long)request);
#endif
#ifdef GA_TOOLS
      char *buffer=(char *)buff;

      int ihandle=(int)*handle;
      int irequest_id;
      *ierr=(fortint)SF_Write(ihandle, *byte_offset, *byte_length, buffer, &irequest_id);
      *request_id=(fortint)irequest_id;
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Sf_write is not available in serial case.\n");
      exit(1);
#endif
   }

/* ******************************************************************************************** *\
   Asynchronous read operation.
   Reads number of bytes to the file identified by handle at location offset.
   Operation is guaranteed to be complete when sf_wait called with request_id argument returns.
\* ******************************************************************************************** */
   void PPIDD_Sf_read(fortint *handle,double *byte_offset,double *byte_length, double *buff,fortint *request_id,fortint *ierr) {
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
      if(MPI_Debug) printf("%5d: In PPIDD_Sf_read  : begin. handle=%d\n",ppidd_sf_rank(),(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      offset=(MPI_Offset)(*byte_offset);
      count=(int)((*byte_length)/sizeof(double));
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_read  : before MPI_File_iread_at. handle=%d,offset=%ld,count=%d\n",ppidd_sf_rank(),(int)mpifhandle,(long)offset,count);
      mpierr=MPI_File_iread_at(mpi_fh,offset,buff,count,MPI_DOUBLE,&request);
      *request_id=(fortint)request;
      *ierr=(fortint)mpierr;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_read  : end. handle=%d,ierr=%d,request_id=%d,request=%ld\n",ppidd_sf_rank(),(int)mpifhandle,(int)*ierr,(int)*request_id,(long)request);
#endif
#ifdef GA_TOOLS
      char *buffer=(char *)buff;

      int ihandle=(int)*handle;
      int irequest_id;
      *ierr=(fortint)SF_Read(ihandle, *byte_offset, *byte_length, buffer, &irequest_id);
      *request_id=(fortint)irequest_id;
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Sf_read is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************************************************ *\
   Blocks the calling process until I/O operation associated with request_id completes.
   Invalidates request_id.
\* ************************************************************************************ */
   void PPIDD_Sf_wait(fortint *request_id,fortint *ierr) {
#ifdef MPI2
      int mpierr;
#ifdef MPIO_USES_MPI_REQUEST
      MPI_Request request=(MPI_Request)(*request_id);
#else
      MPIO_Request request=(MPIO_Request)(*request_id);
#endif
      MPI_Status status;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_wait  : begin. request_id=%d,request=%ld\n",ppidd_sf_rank(),(int)*request_id,(long)request);
#ifdef MPIO_USES_MPI_REQUEST
      mpierr=MPI_Wait( &request, &status );
#else
      mpierr=MPIO_Wait(&request, &status);
#endif
      *ierr=(fortint)mpierr;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_wait  : end. ierr=%d\n",ppidd_sf_rank(),(int)*ierr);
#endif
#ifdef GA_TOOLS
      int irequest_id=(int)*request_id;
      *ierr=(fortint)SF_Wait(&irequest_id);
      *request_id=(fortint)irequest_id;
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Sf_wait is not available in serial case.\n");
      exit(1);
#endif
   }


/* ********************************************************************************** *\
   Blocks the calling process until all of the num I/O operations associated with ids
   specified in list complete. Invalidates (modifies) ids on the list.
\* ********************************************************************************** */
   void PPIDD_Sf_waitall(fortint *list, fortint *num,fortint *ierr) {
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
      *ierr=(fortint)mpierr;
#endif
#ifdef GA_TOOLS
      int inum=(int)*num;
      int *ilist,i;

      ilist=(int *)malloc(inum * sizeof(int));
      for (i=0; i<inum; ++i) {
          ilist[i] = (int) list[i];
      }

      *ierr=(fortint)SF_Waitall(ilist, inum);

      for (i=0; i<inum; ++i) {
          list[i] = (fortint) ilist[i];
      }
      free(ilist);
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Sf_waitall is not available in serial case.\n");
      exit(1);
#endif
   }


/* ************************************************ *\
   Destroys the shared file associated with handle.
   It is a collective operation.
\* ************************************************ */
   void PPIDD_Sf_destroy(fortint *handle,fortint *ierr) {
#ifdef MPI2
      MPI_Fint mpifhandle=(MPI_Fint)*handle;
      MPI_File mpi_fh;
      int mpierr;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_destroy: begin. handle=%d\n",ppidd_sf_rank(),(int)mpifhandle);
      mpi_fh = MPI_File_f2c(mpifhandle);
      mpierr=MPI_File_close( &mpi_fh );
      *ierr=(fortint)mpierr;
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_destroy: end. handle=%d,ierr=%d\n",ppidd_sf_rank(),(int)mpifhandle,(int)*ierr);
#endif
#ifdef GA_TOOLS
      int ihandle=(int)*handle;
      *ierr=(fortint)SF_Destroy(ihandle);
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
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
   void PPIDD_Sf_errmsg(fortint *ierr,char *message
#if defined(FORTCL_NEXT)
	       ,fortintc lx
#endif
#if defined(FORTCL_END)
	       ,fortintc lx
#endif
)  {
#if defined(MPI2) || defined(GA_TOOLS)
#ifdef MPI2
      int eclass, len;
      char estring[MPI_MAX_ERROR_STRING],estring2[MPI_MAX_ERROR_STRING];
#endif
      int perrcode=(int)*ierr;
      int i;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(message);
#endif

      if(MPI_Debug)printf("%5d: In PPIDD_Sf_errmsg: begin. perrcode=%d\n",ppidd_sf_rank(),perrcode);
#ifdef MPI2
      MPI_Error_class(perrcode, &eclass);
      MPI_Error_string(perrcode, estring, &len);
      sprintf(estring2," Error %d: %s", eclass, estring);
      strcpy(message,estring2);
#endif
#ifdef GA_TOOLS
      SF_Errmsg(perrcode, message);
#endif
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_errmsg: middle. message=%s\n",ppidd_sf_rank(),message);
      for(i=strlen(message);i<lxi;i++) message[i]=' ';
      if(MPI_Debug)printf("%5d: In PPIDD_Sf_errmsg: end. perrcode=%d\n",ppidd_sf_rank(),perrcode);
#endif

#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Sf_errmsg is not available in serial case.\n");
      exit(1);
#endif
   }
