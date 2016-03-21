#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

/* ====================================================================== *\
 *                    Standalone Utility Tools for MPI                    *
 *                    ================================                    *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ---------------------------------------------------------------------- *
 * C sorce code of Standalone Utility Tools for MPI. The subroutines can  *
 * be called directly by C code, while the corresponding Fortran wrappers *
 * are unnecessary at the moment.                                         *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       11/03/2009                                                 *
\* ====================================================================== */

#if defined(MPI2) || defined(GA_MPI)


#include "mpi_utils.h"
#include "ppidd_machines.h"  /* needed by fortint in mpiga_type_f2cmpi */

char  mpi_test_err_string[TEST_ERR_STR_LEN];

static int DEBUG_=0;


/* Get the Total node number for specified communicator, and give whether all the nodes are symmetric(every node is the same) */
int NNodes_Total(MPI_Comm comm, int *flag_sym)
{
    int nprocs,rank,nnodes;
    char **nodename;
    int *nprocs_node;
    int length;
    int max_length=256;
    int i,j,skip;
    int sym=1;

    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&rank);

    nodename = (char **)malloc(nprocs * sizeof(char *));
    nprocs_node = (int *)malloc(nprocs * sizeof(int));
    for(i = 0; i < nprocs; i++) {
       nodename[i] = (char *)malloc(max_length * sizeof(char));
       nprocs_node[i] = 0;
    }
#ifdef __bg__
    /* bug 4239 */
    gethostname(nodename[rank],max_length);
    length=strlen(nodename[rank]);
#else
    MPI_Get_processor_name(nodename[rank], &length);
#endif
    if(DEBUG_) fprintf(stdout,"%5d: In NNodes_Total: procname=%s,strlen=%d\n",rank,nodename[rank],length);

    for(i = 0; i < nprocs; i++){
       MPI_Bcast ( nodename[i], max_length, MPI_CHAR, i, comm );
    }

    nnodes=1;
    for(i = 0; i < nprocs; i++){
       skip=0;
       for(j = 0; j < nnodes; j++){
         if(strcmp(nodename[i],nodename[j])==0) skip=1;
       }
       if (skip==0) {nnodes+=1;strcpy(nodename[nnodes-1],nodename[i]);}
       nprocs_node[nnodes-1]+=1;
    }
    /* determine whether all the nodes are symmetric */
    if (nnodes==1) sym=1;
    else {
      for(j = 1; j < nnodes; j++){
        if(nprocs_node[j] != nprocs_node[j-1]) {sym=0;break; }
      }
    }
    *flag_sym=sym;
    if(DEBUG_) fprintf(stdout,"%5d: In NNodes_Total: nnodes=%d, symmetric=%d\n",rank,nnodes,sym);

    for(i = 0; i < nprocs; i++) free(nodename[i]);
    free(nodename);
    free(nprocs_node);

    return(nnodes);
}

/* Get current calling process id */
int ProcID()
{
    int myid;

    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    return(myid);
}


/* Print the error message and exit */
void MPIGA_Error(const char *string, int code)
{
    fprintf(stdout,"%5d: %s %d (%#x).\n", ProcID(), string,code,code);
    fflush(stdout);
    fprintf(stderr,"%5d: %s %d (%#x).\n", ProcID(), string,code,code);

    printf("%5d: In mpi_utils.cpp [MPIGA_Error]: now exiting...\n",ProcID());
#ifdef USE_MPI_ABORT
    MPI_Abort(MPI_COMM_WORLD,code);
#else
    exit(1);
#endif
}


/* Check the returned status; if status is not failing one, then print the error message and exit */
void mpi_test_status(const char *msg_str, int status)
{
    if ( status != MPI_SUCCESS ) {
       int len_err_str, len_msg_str = PPIDD_MIN(MSG_ERR_STR_LEN, strlen(msg_str));
       strncpy(mpi_test_err_string, msg_str, len_msg_str);
       MPI_Error_string(status, mpi_test_err_string + len_msg_str, &len_err_str);
       MPIGA_Error(mpi_test_err_string, status);
    }
}

/* Convert Fortran data type to C MPI_Datatype */
/*
   Input : Fortran Data type ( an int )
   Output: C MPI_Datatype, and size of datatype
===========================================================================
   Fortran data type              (value)    MPI_Datatype
---------------------------------------------------------------------------
   fortint(Integer and Logical)   (0)        MPI_INT,MPI_LONG,MPI_LONG_LONG
   Double Precision               (1)        MPI_DOUBLE
   others                                    not allowed so far, ERROR
===========================================================================
*/
int mpiga_type_f2cmpi(int fdtype, MPI_Datatype *dtype, int *sizeofdtype)
{
    int sizefdtype;
    int sizempidtype;
    MPI_Datatype mpidtype=MPI_CHAR;
    switch(fdtype){
    case 0:
            sizefdtype=sizeof(fortint);
            MPI_Type_size(FORTINT_MPI, &sizempidtype);
            if (sizefdtype==sizempidtype) mpidtype=FORTINT_MPI;
            if (mpidtype==MPI_CHAR) MPIGA_Error("mpiga_type_f2cmpi: can't assign C MPI_Datatype for Fortran Integer ",fdtype);
            break;
    case 1:
            sizefdtype=sizeof(double);
            MPI_Type_size(MPI_DOUBLE, &sizempidtype);
            if (sizefdtype==sizempidtype) mpidtype=MPI_DOUBLE;
            else  MPIGA_Error("mpiga_type_f2cmpi: can't assign C MPI_Datatype for Fortran Double Precision ",fdtype);
            break;
    default:
            MPIGA_Error("mpiga_type_f2cmpi: can't assign C MPI_Datatype for Fortran data type ",fdtype);
            break;
    }
    *dtype=mpidtype;
    *sizeofdtype=sizempidtype;
    return 0;
}


/* Convert C data type to C MPI_Datatype */
/*
   Input : C Data type ( an int )
   Output: C MPI_Datatype, and size of datatype
==========================================================
   MPI_Datatype    C Type       C Constant name   (value)
----------------------------------------------------------
   MPI_INT         int          PPIDD_C_INT       (0)
   MPI_DOUBLE      double       PPIDD_C_DOUBLE    (1)
   MPI_LONG        long         PPIDD_C_LONG      (2)
   MPI_LONG_LONG   long long    PPIDD_C_LONG_LONG (3)
   MPI_FLOAT       float        PPIDD_C_FLOAT     (4)
   others          not allowed so far, ERROR
==========================================================
*/
int mpiga_type_c2cmpi(int cdtype, MPI_Datatype *dtype, int *sizeofdtype)
{
    int scdtype;
    int sizempidtype;
    MPI_Datatype mpidtype=MPI_CHAR;
    switch(cdtype){
    case PPIDD_C_INT:
            scdtype=sizeof(int);
            MPI_Type_size(MPI_INT, &sizempidtype);
            if (scdtype==sizempidtype) mpidtype=MPI_INT;
            break;
    case PPIDD_C_DOUBLE:
            scdtype=sizeof(double);
            MPI_Type_size(MPI_DOUBLE, &sizempidtype);
            if (scdtype==sizempidtype) mpidtype=MPI_DOUBLE;
            break;
    case PPIDD_C_LONG:
            scdtype=sizeof(long);
            MPI_Type_size(MPI_LONG, &sizempidtype);
            if (scdtype==sizempidtype) mpidtype=MPI_LONG;
            break;
    case PPIDD_C_LONG_LONG:
            scdtype=sizeof(long long);
            MPI_Type_size(MPI_LONG_LONG, &sizempidtype);
            if (scdtype==sizempidtype) mpidtype=MPI_LONG_LONG;
            break;
    case PPIDD_C_FLOAT:
            scdtype=sizeof(float);
            MPI_Type_size(MPI_FLOAT, &sizempidtype);
            if (scdtype==sizempidtype) mpidtype=MPI_FLOAT;
            break;
    default:
            MPIGA_Error("mpiga_type_c2cmpi: unknown C data type ",cdtype);
            break;
    }
    if (mpidtype==MPI_CHAR) MPIGA_Error("mpiga_type_c2cmpi: can't convert C data type to C MPI_Datatype ",cdtype);
    *dtype=mpidtype;
    *sizeofdtype=sizempidtype;
    return 0;
}

/* Convert C MPI_Datatype to C data type   */
/*
   Input : C MPI_Datatype
   Output: an int which corresponds to a C Data Type
==========================================================
   MPI_Datatype    C Type       C Constant name   (value)
----------------------------------------------------------
   MPI_INT         int          PPIDD_C_INT       (0)
   MPI_DOUBLE      double       PPIDD_C_DOUBLE    (1)
   MPI_LONG        long         PPIDD_C_LONG      (2)
   MPI_LONG_LONG   long long    PPIDD_C_LONG_LONG (3)
   MPI_FLOAT       float        PPIDD_C_FLOAT     (4)
   others          not allowed so far, ERROR
==========================================================
*/
int mpiga_type_cmpi2c(MPI_Datatype mpidtype, int *cdtype)
{

    if (mpidtype==MPI_INT)            *cdtype=PPIDD_C_INT;
    else if (mpidtype==MPI_LONG)      *cdtype=PPIDD_C_LONG;
    else if (mpidtype==MPI_LONG_LONG) *cdtype=PPIDD_C_LONG_LONG;
    else if (mpidtype==MPI_FLOAT)     *cdtype=PPIDD_C_FLOAT;
    else if (mpidtype==MPI_DOUBLE)    *cdtype=PPIDD_C_DOUBLE;
    else {
       MPIGA_Error("mpiga_type_cmpi2c: unknown C MPI_Datatype ",0);
    }
    return 0;
}


#else

void mpi_utils_dummy () {}

#endif /* end of MPI2 or GA_MPI */
