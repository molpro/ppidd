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

#ifdef HAVE_MPI_H


#include "mpi_utils.h"
#include "ppidd.h"
#include <vector>

char  mpi_test_err_string[TEST_ERR_STR_LEN];

static int DEBUG_=0;

MPI_Datatype dtype_mpi[3];

/* Get the Total node number for specified communicator, and give whether all the nodes are symmetric(every node is the same) */
int NNodes_Total(MPI_Comm comm, int *flag_sym)
{
    int nprocs,rank,nnodes;
    int length;
    constexpr int max_length=256;
    int i,j,skip;
    int sym=1;

    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&rank);

    std::vector<char*> nodename(nprocs);
    nodename[0] = (char*)malloc(nprocs * max_length);
    memset(nodename[0], 0, nprocs * max_length);
    std::vector<int> nprocs_node(nprocs,0);
    for(i = 0; i < nprocs; i++) {
       nodename[i] = nodename[0] + i * max_length;
    }
    MPI_Get_processor_name(nodename[rank], &length);
    if(DEBUG_) fprintf(stdout,"%5d: In NNodes_Total: procname=%s,strlen=%d\n",rank,nodename[rank],length);

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                  nodename[0], max_length, MPI_CHAR, comm);

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

    free(nodename[0]);

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
    exit(1);
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

#else

void mpi_utils_dummy () {}

#endif
