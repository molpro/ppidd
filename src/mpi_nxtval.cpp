#ifdef HAVE_CONFIG_H
#include "ppidd_config.h"
#endif

/* ====================================================================== *\
 *                    MPI Version of Shared Counter                       *
 *                    =============================                       *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ---------------------------------------------------------------------- *
 * C sorce code of MPI Version of Shared Counter. The subroutines can be  *
 * called directly by C code, while the corresponding Fortran wrappers    *
 * (which can be called by Fortran code) are in other files.              *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       02/11/2008                                                 *
\* ====================================================================== */

#ifdef HAVE_MPI_H

#include "ppidd.h"
#include "mpi_utils.h"
#include "mpi_nxtval.h"
#include <cmath>
#include <vector>
int mpiga_cleanup_finalize();

/* twosided helpga global variables */
twosided_helpga_array_t *twosided_helpga_data_struc=NULL, *twosided_helpga_index=NULL;
int *twosided_helpga_map=NULL, *twosided_helpga_proclist=NULL;
int twosided_helpga_num=0;
long twosided_helpga_curmem=0;          /* currently allocated memory of HELPGA for all processes, exclude freed memory */
int HELPGA_LENTOT_SMALL_LIMIT=1024;     /* if the total length of helpga array is less than this limit, then all elements are "allocated" on the first compute process */

/* twosided helpmutex global variables */
twosided_helpmutex_array_t *twosided_helpmutex_index=NULL;            /* twosided helpmutex list */
int twosided_helpmutex_num=0;                                         /* number of twosided helpmutexes */

/* call subroutine PPIDD_Helper_server to change the values of use_helper_server and NPROCS_PER_HELPER */
int use_helper_server=1;         /* helper_server flag: 1 (use); 0 (don't use).                                  */
int NPROCS_PER_HELPER=99999999;  /* how many processes own a helper server.                                      */
int NUMBER_OF_SERVER=1;          /* number of helper servers. Default value will be reset by calling TotalNumber_of_Servers(). */

static int DEBUG_=0;

typedef int64_t dataserver;
#define DATASERVER_MPI MPI_INT64_T

/* Get the total number of server processes */
int TotalNumber_of_Servers()
{
    int numprocs,num;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (use_helper_server) {
      if (numprocs % NPROCS_PER_HELPER ==0 ) num=numprocs/NPROCS_PER_HELPER;
      else num=numprocs/NPROCS_PER_HELPER + 1;
    }
    else {
      num=0;
    }
    if(DEBUG_) fprintf(stdout,"%5d: In TotalNumber_of_Servers: total servers=%d\n",ProcID(),num);
    return num;
}

/* Get the number of work/compute processes */
int NProcs_Work()
{
    int numprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (use_helper_server) {
      if(SR_parallel) return(numprocs-NUMBER_OF_SERVER);
    }
    else {
      if(DEBUG_) fprintf(stdout,"%5d: In NProcs_Work: helper server has not been defined, total procs=%d\n",ProcID(),numprocs);
    }
    return numprocs;
}

/* Get the ID of NXTVAL server: always the last process */
int LastServerID()
{
    int numprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    return(numprocs-1);
}


/* Get the serial number of server process [rank_server]: 0, 1, ..., NUMBER_OF_SERVER-1 */
int SerialNumber_of_Server(int rank_server)
{
    int numprocs,num;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (rank_server == numprocs-1 && numprocs % NPROCS_PER_HELPER!=0) num=(rank_server+1)/NPROCS_PER_HELPER;
    else num=(rank_server+1)/NPROCS_PER_HELPER -1;
    if(DEBUG_) fprintf(stdout,"%5d: In SerialNumber_of_Server: original rank=%d, SerialNumber=%d\n",ProcID(),rank_server,num);
    return num;
}

/* Get the original rank of server process with serial number[0, 1, ..., NUMBER_OF_SERVER-1] */
int RankNumber_of_Server(int rank_serial)
{
    int numprocs,num;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (rank_serial == NUMBER_OF_SERVER-1) num=numprocs-1;
    else num=(rank_serial+1)*NPROCS_PER_HELPER-1;
    if(DEBUG_) fprintf(stdout,"%5d: In RankNumber_of_Server: SerialNumber=%d, original rank=%d\n",ProcID(),rank_serial,num);
    return num;
}

/* Get the corresponding server rank(original) for a compute process(original)*/
int Server_of_Rank(int rank)
{
    int numprocs;
    int server;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    server=(rank/NPROCS_PER_HELPER +1)*NPROCS_PER_HELPER-1;
    if (server > numprocs-1) server = numprocs-1;
    if(DEBUG_) fprintf(stdout,"%5d: In Server_of_Rank: work process rank=%d, server rank=%d\n",ProcID(),rank,server);
    return server;
}

/* Get the number of processes served by server process [rank_server] */
int Nprocs_of_Server(int rank_server)
{
    int numprocs,num;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (rank_server == numprocs-1 && numprocs % NPROCS_PER_HELPER != 0) num=numprocs % NPROCS_PER_HELPER -1;
    else  num=NPROCS_PER_HELPER-1;
    if(DEBUG_) fprintf(stdout,"%5d: In Nprocs_of_Server: server rank=%d, served work processes number=%d\n",ProcID(),rank_server,num);
    return num;
}


/* Get the New ID of old process */
int NewRank_of_OldRank(int rank)
{
    int rank_new;
    int numprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    /* if numprocs==1, then rank_new=rank_old */
    if ( numprocs == 1 ) return rank;
    else {
      if ((rank+1)%NPROCS_PER_HELPER==0 || rank == numprocs-1 ) {
         fprintf(stderr,"%5d: In NewRank_of_OldRank: original process [ %5d ] served as a helper server.\n",ProcID(),rank);
         return 999999;
      }
      else {
         rank_new=rank-rank/NPROCS_PER_HELPER;
         return rank_new;
      }
    }
}

/* Get the original rank of current process in MPIGA_WORK_COMM */
int OldRank_of_NewRank(int rank)
{
    int rank_old;
    int numprocs;
    int numprocs_work=1;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    /* if numprocs==1, then rank_old=rank_new */
    if ( numprocs == 1 ) return rank;
    else {
      numprocs_work=numprocs-NUMBER_OF_SERVER;
      if(rank >  numprocs_work-1) {
         fprintf(stderr,"%5d: In OldRank_of_NewRank: over range rank=%5d, total work processes=%5d.\n",ProcID(),rank,numprocs_work);
         return 999999;
      }
      else {
         rank_old=rank+rank/(NPROCS_PER_HELPER-1);
         return rank_old;
      }
    }
}



/* ===================================================================================== *\
 *                  Create work or computation communicator.                             *
 *  It is a collective operation for old_comm; it will split off one/more processes of   *
 *  the communicator's group to act as data servers, and return a new work communicator. *
 *  The work communicator will be used by the remaining processes for computation.       *
\* ===================================================================================== */
void make_worker_comm( MPI_Comm old_comm, MPI_Comm *worker_comm )
{
    int myid, numprocs, color;

    if(DEBUG_) fprintf(stdout,"%5d: In make_worker_comm begin: use_helper_server=%d\n",ProcID(),use_helper_server);
    MPI_Bcast ( &use_helper_server, 1, MPI_INT, 0, old_comm );
    MPI_Bcast ( &NPROCS_PER_HELPER, 1, MPI_INT, 0, old_comm );
    if (use_helper_server && SR_parallel) {
     /* data servers can be more than a single process */
     MPI_Comm_size(old_comm, &numprocs);
     MPI_Comm_rank(old_comm, &myid);
     if (myid == Server_of_Rank(myid)) color = MPI_UNDEFINED;
     else color = 0;
     MPI_Comm_split( old_comm, color, myid, worker_comm );
    }
    else MPI_Comm_dup( old_comm, worker_comm ); /* duplicate old_com to worker_comm */
    if(DEBUG_) fprintf(stdout,"%5d: In make_worker_comm end: use_helper_server=%d\n",ProcID(),use_helper_server);
}


void DataHelperServer()
{
  int cnt = 0;                 /* actual counter */
  dataserver buf[4];           /* buffer to get values */
  void  *buf_helpga=NULL;      /* help ga which is located in help process */
  void  *buf_ielem =NULL;      /* pointer to helpga[ielem-1]               */
  void  *buf_temp=NULL;        /* temporary buffer which are used to store data for accumulation */
  int32_t *i32buf=NULL,*i32buf_helpga=NULL;
  int64_t *i64buf=NULL,*i64buf_helpga=NULL;
  double *dbuf=NULL,*dbuf_helpga=NULL;
  int  totworkproc;            /* total number of work processes */
  int  Nprocs_server;          /* number of work processes served by current helper server */
  int  type_nxtval = NXTVALFLAG;      /* nxtval operation type            */
  int  type_rma = RMAONEFLAG;         /* one-element RMA operation type   */
  int  type_extra = RMAETRFLAG;       /* extra RMA operation type         */
  int  type_mutex_coll = MUTCOLFLAG;  /* collective mutex operation type  */
  int  type_mutex_lock = MUTLOCFLAG;  /* mutex lock/unlock operation type */
  int  number,inum;
  int  mproc;
  int  opertype;               /* NXTVALFLAG, COLLECFLAG, RMAONEFLAG, RMAETRFLAG, or MUTCOLFLAG */
  int  ndone   = 0;            /* no. finished for this loop */
  int  ntermin=0;
  int  ndone_crt=0, ndone_rls=0, ndone_zero=0;
  int  ndone_mutex_alloc=0, ndone_mutex_free=0;
  int  nodefrom,nodefrom_tmp;
  int  rank_new,rank_next,rank_next_orig;
  int  i,j;
  int  ielem_cdtype,ielem,handle,handle_orig;
  int  nelem=0;
  int  nelem_helpga=0;
  MPI_Request request,reqs[2];
  MPI_Status status,stats[2];
  MPI_Datatype dtype;             /* MPI Datatype */
  int sizeofdtype;                /* Size of MPI Datatype */
  MPIHELPGA helpga=NULL;
  char *p_mutex=NULL;

  if (DEBUG_) printf("%5d: In mpi_nxtval: DataHelperServer begin.\n",ProcID());

  totworkproc=NProcs_Work();
  Nprocs_server=Nprocs_of_Server(ProcID());

  std::vector<int> done_list(totworkproc);
  std::vector<int> done_list_crt(Nprocs_server);
  std::vector<int> done_list_zero(Nprocs_server);
  std::vector<int> done_list_rls(Nprocs_server);
  std::vector<int> done_mutex_alloc(totworkproc);
  std::vector<int> done_mutex_free(totworkproc);

  if (DEBUG_) printf("%5d: In mpi_nxtval: DataHelperServer before while loop. totworkproc=%5d, Nprocs_server=%5d\n",ProcID(),totworkproc,Nprocs_server);

  while (1) {

    /* Wait for input from any process */

    MPI_Recv(buf, 4, DATASERVER_MPI, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    nodefrom = status.MPI_SOURCE;

    opertype = status.MPI_TAG; /* Operartion Flag: NXTVALFLAG (NXTVAL); COLLECFLAG and RMAONEFLAG (twosided_helpga_col and twosided_helpga_one); RMAETRFLAG (twosided_helpga_extra) */
    mproc = (int) buf[0];      /* number of process */


    if (DEBUG_) printf("%5d: DataHelperServer: In the beginning of while loop. from=%d, mproc=%d, opertype=%d\n",ProcID(),nodefrom,mproc,opertype);

    if ( opertype == NXTVALFLAG ) {  /* NXTVAL */
      if (mproc == 0) {  /* terminate helper servers */

      /* Sending process is about to terminate. Send reply and disable
       * sending to him. If all processes have finished return.
       *
       * All processes block on waiting for message
       * from nxtval server before terminating. nxtval only lets
       * everyone go when all have registered termination.
       */
      /* Wait until signals from compute processes, which are served by current server, have finished before terminating it.
       * if Nprocs_server==0, then there is no need to wait */
        done_list[ntermin++] = nodefrom;
        if (ntermin == Nprocs_server || Nprocs_server==0 ) {
          while (ntermin--) {
            nodefrom_tmp = done_list[ntermin];
            MPI_Send(&cnt, 1, MPI_INT,  nodefrom_tmp, type_nxtval, MPI_COMM_WORLD);
          }
          if (DEBUG_) printf("%5d: In mpi_nxtval: DataHelperServer: terminate DataHelperServer.\n",ProcID());
          return;
        }
      }
      else if (mproc > 0) {    /* fetch-and-add INCR */
        MPI_Send(&cnt, 1, MPI_INT,  nodefrom, type_nxtval, MPI_COMM_WORLD);
        cnt += INCR;
      }
      else if (mproc < 0) {   /* release and set cnt=0 */
        /* Wait until all signal from compute processes have finished before releasing it */
        done_list[ndone++] = nodefrom;
        if (ndone == totworkproc) {
          while (ndone--) {
            nodefrom_tmp = done_list[ndone];
            MPI_Send(&cnt, 1, MPI_INT,  nodefrom_tmp, type_nxtval, MPI_COMM_WORLD);
          }
          cnt = 0;
          ndone = 0;
        }
      }
    } /* End of NXTVAL */
    else if ( opertype == COLLECFLAG ) {
       /* HELPGA collective operations: create (mproc=0) , zeroize(mproc>0), destroy (mproc<0) */
      ielem_cdtype=(int)buf[2];     /* COLLECFLAG and mproc=0: C Datatype for the helpga buffer; COLLECFLAG and mproc=others: no use */
      handle=(int)buf[3];           /* adjuested sequence number of helpga */
      handle_orig=handle-TWOSIDED_HELPGA_OFFSET; /* original sequence number of helpga */

      if (mproc == 0) { /* create a help GA */

        /* All work processes tell server to create an array with buf[1] elements.
         * Servers wait until it has received all requests from all work processes.
         * After received all requests, helpga is created in server, and server send
         * a SUCESS flag to all work processes.
         *
         * All work processes block on waiting for message from nxtval server before return.
         * If all work processes have received SUCESS flag, then return.
         */
        /* Wait until signals from compute processes, which are served by current server, have finished before creating it */
        done_list_crt[ndone_crt++] = nodefrom;
        if (ndone_crt == Nprocs_server) {
          nelem_helpga=(int)buf[1];
          dtype=dtype_mpi[ielem_cdtype];
          MPI_Type_size( dtype, &sizeofdtype );

          /* create  a new structure */
          helpga =(MPIHELPGA)malloc( sizeof(struct STRUC_MPIHELPGA) );
          if (!helpga) MPIGA_Error("DataHelperServer ERROR: failed to allocate memory for STRUC_MPIHELPGA ",0);

          helpga->name=NULL;
          helpga->dtype=dtype;
          helpga->nele=nelem_helpga;        /* total number of elements on current server */
          helpga->len=NULL;
          helpga->len_help=NULL;

/*          if(nelem_helpga!=0) helpga->ptr_buf=malloc(nelem_helpga*sizeofdtype); */
          if(nelem_helpga!=0) MPI_Alloc_mem(nelem_helpga*sizeofdtype, MPI_INFO_NULL, &helpga->ptr_buf);
          else helpga->ptr_buf=NULL;

          twosided_helpga_index[handle_orig].actv=1;
          twosided_helpga_index[handle_orig].ptr=helpga;
          twosided_helpga_num++;
          twosided_helpga_curmem += nelem_helpga*sizeofdtype;

          while (ndone_crt--) {
            nodefrom_tmp = done_list_crt[ndone_crt];
            MPI_Send(NULL, 0, MPI_BYTE,  nodefrom_tmp, WAKEUPTAG, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
          }
          ndone_crt = 0;
          helpga=NULL;
          if (DEBUG_) printf("%5d: DataHelperServer: helpga_create end. handle=%d, nelem_helpga=%d, twosided_helpga_num=%d\n",ProcID(),handle,nelem_helpga,twosided_helpga_num);
        }
      } /* end of creating help GA */
      else if (mproc > 0) { /* zeroize help GA */
        /* Wait until signals from compute processes, which are served by current server, have finished before zeroizing it */
        done_list_zero[ndone_zero++] = nodefrom;
        if (ndone_zero == Nprocs_server) {

          while (ndone_zero--) {
            nodefrom_tmp = done_list_zero[ndone_zero];
            MPI_Send(NULL, 0, MPI_BYTE,  nodefrom_tmp, WAKEUPTAG, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
          }

          helpga=twosided_helpga_index[handle_orig].ptr;
          nelem_helpga=helpga->nele;
          dtype=helpga->dtype;
          if (nelem_helpga!=0) {
            if (dtype==MPI_INT32_T) {
               i32buf=(int32_t *)helpga->ptr_buf;
               for (i=0;i<nelem_helpga;i++) i32buf[i]=(int32_t)0;
               i32buf=NULL;
            }
            else if (dtype==MPI_INT64_T) {
               i64buf=(int64_t *)helpga->ptr_buf;
               for (i=0;i<nelem_helpga;i++) i64buf[i]=(int64_t)0;
               i64buf=NULL;
            }
            else if (dtype==MPI_DOUBLE) {
               dbuf=(double *)helpga->ptr_buf;
               for (i=0;i<nelem_helpga;i++) dbuf[i]=0.0e0;
               dbuf=NULL;
            }
            else {
               MPIGA_Error("DataHelperServer [helpga_zero]: wrong MPI_Datatype ",0);
            }
          }
          helpga=NULL;
          ndone_zero = 0;
          if (DEBUG_) printf("%5d: DataHelperServer: helpga_zero end. handle=%d, nelem_helpga=%d, twosided_helpga_num=%d\n",ProcID(),handle,nelem_helpga,twosided_helpga_num);
        }
      } /* end of zeroizing HELP GA */
      else if (mproc < 0) { /* release/destroy HELP GA */
        /* Wait until signals from compute processes, which are served by current server, have finished before releasing it */
        done_list_rls[ndone_rls++] = nodefrom;
        if (ndone_rls == Nprocs_server) {
          while (ndone_rls--) {
            nodefrom_tmp = done_list_rls[ndone_rls];
            MPI_Send(NULL, 0, MPI_BYTE,  nodefrom_tmp, WAKEUPTAG, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
          }
          twosided_helpga_release_orig(handle_orig);
          ndone_rls = 0;
          if (DEBUG_) printf("%5d: DataHelperServer: helpga_release end. handle=%d, nelem_helpga=%d, twosided_helpga_num=%d\n",ProcID(),handle,nelem_helpga,twosided_helpga_num);
        }
      } /* release/destroy HELP GA */
    } /* End of collective HELPGA */
    else if ( opertype == RMAONEFLAG ) {
       /* HELPGA RMA operations: get (mproc=0) , fetch-and-add(mproc>0), put (mproc<0) */
      ielem =(int)buf[2];           /* sequence number of element in helpga: 1,2,...,n      */
      handle=(int)buf[3];           /* adjuested sequence number of helpga */
      handle_orig=handle-TWOSIDED_HELPGA_OFFSET; /* original sequence number of helpga */

      helpga=twosided_helpga_index[handle_orig].ptr;
      dtype=helpga->dtype;
      if (dtype==MPI_INT32_T) {
         i32buf = (int32_t *)helpga->ptr_buf;
         buf[0] = (dataserver)i32buf[ielem-1];
     /*  if (mproc == 0) work process gets a value from help GA, ie server sends an element value to a work process */
         if (mproc > 0) i32buf[ielem-1] += (int32_t)buf[1]; /* fetch-and-add INCR to help GA */
	 if (mproc < 0) i32buf[ielem-1] = (int32_t)buf[1]; /* put a value to help GA */
      }
      else if (dtype==MPI_INT64_T) {
         i64buf = (int64_t *)helpga->ptr_buf;
         buf[0] = (dataserver)i64buf[ielem-1];
     /*  if (mproc == 0) work process gets a value from help GA, ie server sends an element value to a work process */
         if (mproc > 0) i64buf[ielem-1] += (int64_t)buf[1]; /* fetch-and-add INCR to help GA */
	 if (mproc < 0) i64buf[ielem-1] = (int64_t)buf[1]; /* put a value to help GA */
      }
      MPI_Send(buf, 1, DATASERVER_MPI,  nodefrom, type_rma, MPI_COMM_WORLD);
      if (DEBUG_) printf("%5d: DataHelperServer: helpga_one end. handle=%d, ielem_helpga=%d, twosided_helpga_num=%d, returnval=%ld\n",ProcID(),handle,ielem,twosided_helpga_num,(long)buf[0]);
      helpga=NULL;
    } /* End of HELPGA RMA operations */
    else if ( opertype == RMAETRFLAG ) {
       /* HELPGA extra RMA operations (RMAETRFLAG): get(mproc=0), put(mproc<0), accumulate(mproc>0) */
      nelem_helpga=(int)buf[1];     /* mproc<=0: number of elements to be gotten/put/accumulated  */
      ielem =(int)buf[2];           /* sequence number of element in helpga: 1,2,...,n      */
      handle=(int)buf[3];           /* adjuested sequence number of helpga */
      handle_orig=handle-TWOSIDED_HELPGA_OFFSET; /* original sequence number of helpga */

      helpga=twosided_helpga_index[handle_orig].ptr;
      dtype=helpga->dtype;
      MPI_Type_size( dtype, &sizeofdtype );

      if (mproc == 0) { /* work process gets a set of values from help GA, ie server sends a set of values to a work process */
        buf_helpga=helpga->ptr_buf;
        buf_ielem = (void *)( (char *)buf_helpga + (ielem-1) * sizeofdtype );
        MPI_Isend(buf_ielem, nelem_helpga, dtype,  nodefrom, type_extra, MPI_COMM_WORLD, &request);  /*send n elements to work process*/
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        if (DEBUG_) MPI_Send(NULL, 0, MPI_BYTE,  nodefrom, WAKEUPTAG, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
      }
      else if (mproc < 0) { /* put a set of values to help GA */
        /*receive n elements from work process*/
        /* method 1: */
        buf_helpga=helpga->ptr_buf;
        buf_ielem = (void *)( (char *)buf_helpga + (ielem-1) * sizeofdtype );
        MPI_Irecv(buf_ielem, nelem_helpga, dtype,  nodefrom, type_extra, MPI_COMM_WORLD, &reqs[0]);
        /* method 2:
        if (!buf_temp) buf_temp=malloc(nelem_helpga*sizeofdtype);
        MPI_Irecv(buf_temp, nelem_helpga, dtype,  nodefrom, type_extra, MPI_COMM_WORLD, &reqs[0]);
        */
        /* notification is necessary to ensure global data structure is updated consistently, otherwise some jobs fail on multinodes */
        MPI_Isend(NULL, 0, MPI_BYTE,  nodefrom, WAKEUPTAG, MPI_COMM_WORLD, &reqs[1]);  /* SUCCESS: notify compute process */
        MPI_Waitall(2, reqs, stats);
        if (DEBUG_) {
          MPI_Get_count(&stats[0], dtype, &nelem);
          if ( nelem != nelem_helpga) MPIGA_Error("DataHelperServer: for put operation. nelem != nelem_helpga", nelem);
        }

        /* method 2:
        if (dtype==MPI_INT) {
           ibuf=(int *)buf_temp;
           ibuf_helpga=(int *)helpga->ptr_buf;
           for (i=0;i<nelem_helpga;i++) ibuf_helpga[ielem-1+i]=ibuf[i];
           ibuf=NULL;
           ibuf_helpga=NULL;
        }
        else if (dtype==MPI_DOUBLE) {
           dbuf=(double *)buf_temp;
           dbuf_helpga=(double *)helpga->ptr_buf;
           for (i=0;i<nelem_helpga;i++) dbuf_helpga[ielem-1+i]=dbuf[i];
           dbuf=NULL;
           dbuf_helpga=NULL;
        }
        else {
           MPIGA_Error("DataHelperServer [twosided_helpga_extra put]: wrong MPI_Datatype ",0);
        }
        if (buf_temp != NULL) { free (buf_temp); buf_temp=NULL; } */
      }
      else { /* mproc>0, accumulate a set of values to help GA (First receive data, then add it to help GA ) */
       /*receive n elements from work process*/
        if (!buf_temp) buf_temp=malloc(nelem_helpga*sizeofdtype);
        MPI_Irecv(buf_temp, nelem_helpga, dtype,  nodefrom, type_extra, MPI_COMM_WORLD, &reqs[0]);
        /* data servers send signals to ensure accumulate is an atomic operation */
        MPI_Isend(NULL, 0, MPI_BYTE,  nodefrom, WAKEUPTAG, MPI_COMM_WORLD, &reqs[1]);  /* SUCCESS: notify compute process */
        MPI_Waitall(2, reqs, stats);
        if (DEBUG_) {
          MPI_Get_count(&stats[0], dtype, &nelem);
          if ( nelem != nelem_helpga) MPIGA_Error("DataHelperServer: for accumulation. nelem != nelem_helpga", nelem);
        }
        if (dtype==MPI_INT32_T) {
           i32buf=(int32_t *)buf_temp;
           i32buf_helpga=(int32_t *)helpga->ptr_buf;
           for (i=0;i<nelem_helpga;i++) i32buf_helpga[ielem-1+i]+=i32buf[i];
           i32buf=NULL;
           i32buf_helpga=NULL;
        }
        else if (dtype==MPI_INT64_T) {
           i64buf=(int64_t *)buf_temp;
           i64buf_helpga=(int64_t *)helpga->ptr_buf;
           for (i=0;i<nelem_helpga;i++) i64buf_helpga[ielem-1+i]+=i64buf[i];
           i64buf=NULL;
           i64buf_helpga=NULL;
        }
        else if (dtype==MPI_DOUBLE) {
           dbuf=(double *)buf_temp;
           dbuf_helpga=(double *)helpga->ptr_buf;
           for (i=0;i<nelem_helpga;i++) dbuf_helpga[ielem-1+i]+=dbuf[i];
           dbuf=NULL;
           dbuf_helpga=NULL;
        }
        else {
           MPIGA_Error("DataHelperServer [twosided_helpga_extra acc]: wrong MPI_Datatype ",0);
        }
        if (buf_temp != NULL) { free (buf_temp); buf_temp=NULL; }
      }
      helpga=NULL;
    } /* End of HELPGA extra RMA operations */
    else if ( opertype == MUTCOLFLAG ) {
       /* twosided helpmutex allocate and free operations: allocate(mproc>0), free(mproc<0),no-use(mproc=0) */
      inum=(int)buf[1];          /* MUTCOLFLAG(inum): first sequential number of mutexes to be allocated (mproc > 0) or freed (mproc < 0)  */
      number=(int)buf[2];        /* MUTCOLFLAG(number): number of mutexes to be allocated (mproc > 0) or freed (mproc < 0) */
      if (mproc > 0) { /* allocate twosided helpmutexes  */
        /* Wait until all signal from compute processes have finished before allocating a set of mutexes */
        done_mutex_alloc[ndone_mutex_alloc++] = nodefrom;
        if (ndone_mutex_alloc == totworkproc) {
          while (ndone_mutex_alloc--) {
            nodefrom_tmp = done_mutex_alloc[ndone_mutex_alloc];
            MPI_Send(NULL, 0, MPI_BYTE,  nodefrom_tmp, type_mutex_coll, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
          }
          /* allocate memory for mutexes, and zeroize them. One mutex is controlled by one int array with length of totworkproc. */
          for(i=inum;i<inum+number; i++) {
            twosided_helpmutex_index[i].actv = 1;
            twosided_helpmutex_index[i].mutex = (char *)malloc(totworkproc*sizeof(char));
            for(j=0;j<totworkproc; j++) (twosided_helpmutex_index[i]).mutex[j]=0;
          }
          if (inum>=MAX_TWOSIDED_HELPMUTEX_GA) twosided_helpmutex_num=number;
          ndone_mutex_alloc = 0;
          if (DEBUG_) printf("%5d: DataHelperServer: twosided_helpmutex_collect end. %d mutexes have been allocated. inum=%d.\n",ProcID(),number,inum);
        }
      } /* allocate twosided helpmutexes */
      else if (mproc < 0) { /* free twosided helpmutexes  */
        /* Wait until all signal from compute processes have freed before allocating a set of mutexes */
        done_mutex_free[ndone_mutex_free++] = nodefrom;
        if (ndone_mutex_free == totworkproc) {
          while (ndone_mutex_free--) {
            nodefrom_tmp = done_mutex_free[ndone_mutex_free];
            MPI_Send(NULL, 0, MPI_BYTE,  nodefrom_tmp, type_mutex_coll, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
          }
          /* free the memory mutexes have occupied */
          for(i=inum;i<inum+number; i++) {
            twosided_helpmutex_index[i].actv = 0;
            if ( twosided_helpmutex_index[i].mutex != NULL ) { free (twosided_helpmutex_index[i].mutex); twosided_helpmutex_index[i].mutex=NULL; }
          }
          if (inum>=MAX_TWOSIDED_HELPMUTEX_GA) twosided_helpmutex_num=0;
          ndone_mutex_free = 0;
          if (DEBUG_) printf("%5d: DataHelperServer: twosided_helpmutex_collect end. %d mutexes have been freed. inum=%d.\n",ProcID(),number,inum);
        }
      } /* free twosided helpmutexes */
      else { /* mproc==0 */
        MPIGA_Error("DataHelperServer: twosided_helpmutex_collect with wrong mproc ", mproc);
      }
    } /* End of twosided helpmutex allocate and free operations */
    else if ( opertype == MUTLOCFLAG ) {
       /* twosided helpmutex lock and unlock operations: lock(mproc>0), unlock(mproc<0),no-use(mproc=0) */
      inum=(int)buf[1]; /* original sequence number of twosided helpmutex */
      rank_new=NewRank_of_OldRank(nodefrom);
      p_mutex=twosided_helpmutex_index[inum].mutex;
      if (mproc > 0) { /* lock the twosided helpmutex numbered inum  */
        /* add self to the waiting list */
        if ( p_mutex[rank_new] != 0 ) MPIGA_Error("DataHelperServer: twosided_helpmutex_lock. attempt to lock a locked mutex by this process. nodefrom ", nodefrom);
        /* check if other compute processes are waiting for the lock */
        for (i=0; i<totworkproc && p_mutex[i]==0; i++);
        /* check if rank_new is the only waiting compute process. If so, then send SUCCESS flag to compute process */
        if ( i >= totworkproc ) {
          MPI_Send(NULL, 0, MPI_BYTE,  nodefrom, type_mutex_lock, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
          if (DEBUG_) printf("%5d: DataHelperServer: twosided_helpmutex_lock end. mutex %d has been locked by process %d.\n",ProcID(),inum,nodefrom);
        }
        /* add self to the waiting list */
        p_mutex[rank_new]=1;
      } /* lock the twosided helpmutex numbered inum */
      else if (mproc < 0) { /* unlock the twosided helpmutex numbered inum */
        /* remove self from the waiting list, and send SUCCESS flag to compute process */
        if (p_mutex[rank_new]==1) p_mutex[rank_new]=0;
        else MPIGA_Error("DataHelperServer: twosided_helpmutex_lock. attempt to unlock a mutex which is not locked by this process. nodefrom ", nodefrom);
        MPI_Send(NULL, 0, MPI_BYTE,  nodefrom, type_mutex_lock, MPI_COMM_WORLD);  /* SUCCESS: notify compute process */
        if (DEBUG_) printf("%5d: DataHelperServer: twosided_helpmutex_lock end. mutex %d has been unlocked by process %d.\n",ProcID(),inum,nodefrom);
        /* check if there is a compute process waiting for the lock */
        for (i=0; i<totworkproc && p_mutex[i]==0; i++);
        if (i<totworkproc) {
        /* find the next waiting compute process, and then send SUCCESS flag to it */
          rank_next=rank_new;
          /* find the next rank who is waiting for the lock */
          while ( rank_next < totworkproc && p_mutex[rank_next] == 0 ) rank_next++;
          if ( rank_next == totworkproc ) {
             rank_next=0;
             while ( rank_next < rank_new && p_mutex[rank_next] == 0 ) rank_next++;
          }
          /* notify the next rank who is waiting for the lock */
          rank_next_orig=OldRank_of_NewRank(rank_next);
          MPI_Send(NULL, 0, MPI_BYTE,  rank_next_orig, type_mutex_lock, MPI_COMM_WORLD);  /* SUCCESS: notify the next rank who is waiting for the lock */
          if (DEBUG_) printf("%5d: DataHelperServer: twosided_helpmutex_lock end. mutex %d has been locked by process %d.\n",ProcID(),inum,rank_next_orig);
        }
      } /* unlock the twosided helpmutex numbered inum */
      else { /* mproc==0 */
        MPIGA_Error("DataHelperServer: twosided_helpmutex_lock with wrong mproc ", mproc);
      }
      p_mutex=NULL;
    } /* End of twosided helpmutex lock and unlock operations */
    else { /* other unknown opertype */
      MPIGA_Error("DataHelperServer: unknown opertype ", opertype);
    }
    if (DEBUG_) printf("%5d: DataHelperServer: In the end of while loop. from=%d, mproc=%d, opertype=%d\n",ProcID(),nodefrom,mproc,opertype);
  } /* end of while loop */
}


int NXTVAL(int mproc)
/*
  Get next value of shared counter.

  mproc > 0 ... fetch-and-add operation, returns requested value
  mproc < 0 ... server blocks until abs(mproc) processes are queued
                and returns junk, release cnt
  mproc = 0 ... indicates to server that I am about to terminate

*/
{
  int  type = NXTVALFLAG;
  int  cnt;
  dataserver buf;
  int  local=0;
  MPI_Status status;
  int  server;                    /* id of server process */
  int  myid;                      /* id of current process */


  if (SR_parallel) {

     buf = mproc;

     if (DEBUG_) {
       printf("%5d: NXTVAL: mproc=%d\n",ProcID(), mproc);
       fflush(stdout);
     }

    if (use_helper_server) {
       MPI_Comm_rank(MPI_COMM_WORLD, &myid);
       if (mproc == 0) server=Server_of_Rank(myid); /* terminate all helper servers */
       else server = LastServerID();                 /* last server process */

       MPI_Send(&buf, 1, DATASERVER_MPI,  server, type, MPI_COMM_WORLD);
       MPI_Recv(&cnt, 1, MPI_INT,  server, type, MPI_COMM_WORLD, &status);
       if (mproc == 0 && myid == 0 ) { /* rank(0) sends signal to terminate last server process if the last server process doesn't serve any compute process */
         if ( Nprocs_of_Server(LastServerID()) == 0 ) {
           server = LastServerID();                 /* last server process */
	   buf = mproc;
           MPI_Send(&buf, 1, DATASERVER_MPI,  server, type, MPI_COMM_WORLD);
           MPI_Recv(&cnt, 1, MPI_INT,  server, type, MPI_COMM_WORLD, &status);
         }
       }
       return cnt;
    }
  }
  else {
     /* Not running in parallel ... just do a simulation */
     static int count = 0;
     if (mproc == 1)
       local = count++;
     else if (mproc == -1) {
       count = 0;
       local = 0;
    }
    else
      MPIGA_Error("NXTVAL: sequential version with silly mproc ", mproc);
  }

  return local;
}

/* Collective operations(zeroize, destroy) of helpga */
int twosided_helpga_col(int mproc, int handle)
/*
  Operations for helpga:
________________________________________________________________________________________________________
|                |     COLLECFLAG       |       RMAONEFLAG             |        RMAETRFLAG             |
|  mproc = 0 ... |   create helpga      | get one element from helpga  | get elements from helpga      |
|  mproc > 0 ... |   zeroize helpga     | fetch-and-add to helpga      | accumulate helpga elements    |
|  mproc < 0 ... |   destroy helpga     | put one value to helpga      | put values to helpga elements |
--------------------------------------------------------------------------------------------------------
*/
{
  int  type = COLLECFLAG;
  dataserver buf[4];
  int  local=0;
  int  handle_orig;
  MPI_Datatype dtype;             /* MPI Datatype */
  int  myid;                      /* id of server compute process */
  int  server;                    /* id of server process         */
  MPIHELPGA helpga=NULL;

  if (DEBUG_) printf("%5d: twosided_helpga_col: begin. type=%d, mproc=%d, handle=%d\n",ProcID(),type,mproc,handle);

  if ( MPIGA_WORK_COMM != MPI_COMM_NULL) MPI_Barrier(MPIGA_WORK_COMM);

  if( mproc==0 ) { /* can't be create operation */
    fprintf(stderr," twosided_helpga_col ERROR:  not zeroize/release operation.\n");
    return 1;
  }
  else { /* zeroize/release operation */

    handle_orig=twosided_helpga_handle_orig(handle);
    helpga=twosided_helpga_index[handle_orig].ptr;

  } /* other operations (zeroization and release) for helpga */

  if (DEBUG_) {
    printf("%5d: twosided_helpga_col: mid. type=%d, mproc=%d, handle=%d\n",ProcID(),type,mproc,handle);
    fflush(stdout);
  }

  if (SR_parallel) {
     buf[0] = (dataserver)mproc; /* COLLECFLAG: create(=0), zeroize(>0), destroy(<0); RMAONEFLAG: get(=0), fetch-and-add(>0), put(<0) */
     buf[1] = (dataserver)1;     /* COLLECFLAG (mproc=0: number of elements; others: no use). RMAONEFLAG: value to be put(mproc<0); increment value(mproc>0); no use(mproc=0)  */
     buf[2] = (dataserver)1;     /* COLLECFLAG (mproc=0: C data type; others: no use).        RMAONEFLAG: sequence number of element (1,2,...,n) */
     buf[3] = (dataserver)handle;      /* sequence number of helpga */

     if (use_helper_server) {
       MPI_Comm_rank(MPI_COMM_WORLD, &myid);
       server=Server_of_Rank(myid);

       MPI_Send(buf, 4, DATASERVER_MPI, server, type, MPI_COMM_WORLD);
       MPI_Recv(NULL, 0, MPI_BYTE, server, WAKEUPTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* SUCCESS: receive notification from server */
     }
  }
  else {
     /* sequential case: Not running in parallel ... just do a simulation */
    int32_t *i32buf;
    int64_t *i64buf;
    double *dbuf;
    int j;
    int nelem_localga=0;
    dtype=helpga->dtype;
    nelem_localga=helpga->nele;
    if (mproc == 1) {            /* zeroize local helpga */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)helpga->ptr_buf;
         for (j=0;j<nelem_localga;j++) i32buf[j]=(int32_t)0;
         i32buf=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)helpga->ptr_buf;
         for (j=0;j<nelem_localga;j++) i64buf[j]=(int64_t)0;
         i64buf=NULL;
      }
      else if (dtype==MPI_DOUBLE) {
         dbuf=(double *)helpga->ptr_buf;
         for (j=0;j<nelem_localga;j++) dbuf[j]=0.0e0;
         dbuf=NULL;
      }
      else {
         MPIGA_Error("twosided_helpga_col [helpga_zero]: wrong MPI_Datatype ",0);
      }
      local = 0;
    }
    else if (mproc == -1) {           /* release local helpga */
      local = 0;
    }
    else {
      MPIGA_Error("twosided_helpga_col: sequential version with silly mproc ", mproc);
    }
  } /* sequential case */

  if (DEBUG_) printf("%5d: twosided_helpga_col: before end. type=%d, mproc=%d, handle=%d\n",ProcID(),type,mproc,handle);

  if( mproc<0 ) { /* release helpga */
    twosided_helpga_release_orig(handle_orig);
  }

  helpga=NULL;

  if ( MPIGA_WORK_COMM != MPI_COMM_NULL) MPI_Barrier(MPIGA_WORK_COMM);

  if (DEBUG_) printf("%5d: twosided_helpga_col: end. type=%d, mproc=%d, handle=%d\n",ProcID(),type,mproc,handle);

  return local;
}


int twosided_helpga_create_irreg(int mproc, int *lenin, int nchunk, int *handle, char *name, int dtype)
/*
  Operations for helpga:
________________________________________________________________________________________________________
|                |     COLLECFLAG       |       RMAONEFLAG             |        RMAETRFLAG             |
|  mproc = 0 ... |   create helpga      | get one element from helpga  | get elements from helpga      |
|  mproc > 0 ... |   zeroize helpga     | fetch-and-add to helpga      | accumulate helpga elements    |
|  mproc < 0 ... |   destroy helpga     | put one value to helpga      | put values to helpga elements |
--------------------------------------------------------------------------------------------------------
*/
{
  int  type = COLLECFLAG;
  dataserver buf[4];
  int  local=0;
  int  handle_orig=0;
  int sizeofdtype;                /* Size of MPI Datatype */
  int  myid;                      /* id of compute process        */
  int  server;                    /* id of server process         */
  int  lentot;                    /* total number of elements     */
  long sizetot;
  int  i,ii;
  int  size,rank_pre,rank_serial;
  char *helpganame=NULL;
  int  *len=NULL;                 /* size on each compute process */
  int  *len_help=NULL;            /* size on each helper process  */
  MPIHELPGA helpga=NULL;

  if (DEBUG_) printf("%5d: twosided_helpga_create_irreg: begin. type=%d, mproc=%d\n",ProcID(),type,mproc);

  if ( MPIGA_WORK_COMM != MPI_COMM_NULL) MPI_Barrier(MPIGA_WORK_COMM);

  if( mproc!=0 ) { /*check if create operation */
    fprintf(stderr," twosided_helpga_create_irreg ERROR:  not create operation.\n");
    return 1;
  }
  else { /*create helpga */
    twosided_helpga_num++;
   /* Check to ensure the maximum number of arrays hasn't been reached.*/
    if( twosided_helpga_num > MAX_TWOSIDED_HELPGA_ARRAYS ) {
      fprintf(stderr," twosided_helpga_create_irreg ERROR:  %d over the maximum number of mpi helpga [%i].\n",twosided_helpga_num,MAX_TWOSIDED_HELPGA_ARRAYS);
      return 1;
    }
    for(i=0;i<MAX_TWOSIDED_HELPGA_ARRAYS;i++){
       if ( twosided_helpga_index[i].actv == 0 ) {
          handle_orig=i;  /* original sequence number of helpga */
          break;
       }
    }

    /* create  a new structure */
    helpga =(MPIHELPGA)malloc( sizeof(struct STRUC_MPIHELPGA) );
    if (!helpga) return 1;

/* calculate the number of elements for each compute process */
    size=NProcs_Work();

    if (nchunk<0 || nchunk>size) {
       printf("ERROR in mpiga_create_irreg : nchunk(%5d)larger than process number(%5d).\n", nchunk,size);
       return 1;
    }

    len=(int *)malloc( size*sizeof(int) );
    lentot=0;
    for(i=0;i<nchunk;i++) {len[i]=lenin[i];lentot=lentot+lenin[i];}
    for(i=nchunk;i<size;i++) len[i]=0;


/* calculate the number of elements for each helper server */
    len_help=(int *)malloc(NUMBER_OF_SERVER*sizeof(int));
    for(i=0; i<NUMBER_OF_SERVER; i++) len_help[i]=0;
    rank_pre=Server_of_Rank(0);
    ii=0;
    for(i=0;i<size;i++) {
       if (Server_of_Rank(OldRank_of_NewRank(i))==rank_pre) len_help[ii]+=len[i];
       else {
         ii++;
         len_help[ii]+=len[i];
         rank_pre=Server_of_Rank(OldRank_of_NewRank(i));
       }
    }


    strcpy(helpganame=(char *)malloc(strlen(name)+1),name);

    MPI_Datatype mpidtype=dtype_mpi[dtype];

    helpga->name=helpganame;
    helpga->ptr_buf=NULL;
    helpga->dtype=mpidtype;
    helpga->nele=lentot;        /* total number of elements */
    helpga->len=len;            /* number of elements for each compute process */
    helpga->len_help=len_help;  /* number of elements for each helper server */

    MPI_Type_size( mpidtype, &sizeofdtype );
    sizetot=(long)sizeofdtype*(long)lentot;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    rank_serial=NewRank_of_OldRank(myid);

    twosided_helpga_index[handle_orig].size=(long)sizeofdtype*(long)len[rank_serial];
    twosided_helpga_index[handle_orig].actv=1;
    twosided_helpga_index[handle_orig].ptr=helpga;

    twosided_helpga_curmem += sizetot;

    *handle=handle_orig+TWOSIDED_HELPGA_OFFSET;  /* adjusted sequence number of helpga */

    len=NULL;
    len_help=NULL;
    helpganame=NULL;

  }

  if (DEBUG_) printf("%5d: twosided_helpga_create_irreg: mid. type=%d, mproc=%d,handle=%d\n",ProcID(),type,mproc,*handle);

  if (SR_parallel) {
     buf[0] = (dataserver)mproc; /* COLLECFLAG: create(=0), zeroize(>0), destroy(<0); RMAONEFLAG: get(=0), fetch-and-add(>0), put(<0) */
     buf[1] = (dataserver)lentot;       /* COLLECFLAG and mproc=0: number of elements; RMAONEFLAG: value to be put(mproc<0), increment value(mproc>0); others: no use */
     buf[2] = (dataserver)dtype;        /* COLLECFLAG and mproc=0: MPI_Datatype; RMAONEFLAG: sequence number of element (1,2,...,n) */
     buf[3] = (dataserver)*handle;      /* sequence number of helpga */

     if (use_helper_server) {
       MPI_Comm_rank(MPI_COMM_WORLD, &myid);
       server=Server_of_Rank(myid);
       ii=SerialNumber_of_Server(server);
       buf[1] = (dataserver) (helpga->len_help[ii]); /* COLLECFLAG and mproc=0: number of elements; RMAONEFLAG: value to be put(mproc<0), increment value(mproc>0); others: no use */

       MPI_Send(buf, 4, DATASERVER_MPI,  server, type, MPI_COMM_WORLD);
       MPI_Recv(NULL, 0, MPI_BYTE, server, WAKEUPTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* SUCCESS: receive notification from server */
     }
  }
  else {
     /* sequential case: Not running in parallel ... just do a simulation */
    int nelem_localga=0;
    nelem_localga=helpga->nele;
    if (mproc == 0) {                 /* create a local helpga  */
/*      if(nelem_localga!=0) helpga->ptr_buf=malloc(nelem_localga*sizeofdtype); */
      if(nelem_localga!=0) MPI_Alloc_mem(nelem_localga*sizeofdtype, MPI_INFO_NULL, &helpga->ptr_buf);
      local = 0;
    }
    else {
      MPIGA_Error("twosided_helpga_create_irreg: sequential version with silly mproc ", mproc);
    }
  } /* sequential case */

  helpga=NULL;

  if ( MPIGA_WORK_COMM != MPI_COMM_NULL) MPI_Barrier(MPIGA_WORK_COMM);

  if (DEBUG_) printf("%5d: twosided_helpga_create_irreg: end. type=%d, mproc=%d,handle=%d\n",ProcID(),type,mproc,*handle);

  return local;
}


/* create a normal helpga */
/* here  mproc = 0 ;  lentot: total number of elements     */
int twosided_helpga_create(int mproc, int lentot, int *handle, char *name, int dtype)
/*
  Operations for helpga:
________________________________________________________________________________________________________
|                |     COLLECFLAG       |       RMAONEFLAG             |        RMAETRFLAG             |
|  mproc = 0 ... |   create helpga      | get one element from helpga  | get elements from helpga      |
|  mproc > 0 ... |   zeroize helpga     | fetch-and-add to helpga      | accumulate helpga elements    |
|  mproc < 0 ... |   destroy helpga     | put one value to helpga      | put values to helpga elements |
--------------------------------------------------------------------------------------------------------
*/
{
    int sizeoflen;
    std::vector<int> len; /* size on each compute process */

    if (DEBUG_) printf("%5d: twosided_helpga_create: begin. mproc=%d, lentot=%d\n",ProcID(),mproc,lentot);

/* calculate the number of elements for each compute process */
    int size=NProcs_Work();

    /* len: non-zero array that stores the number of elements on i-th compute process; sizeoflen: size of array len */
    /* if array is very small, then all elements are allocated on the first compute process */
    if (lentot <= HELPGA_LENTOT_SMALL_LIMIT) {
      sizeoflen=1;
      len.resize(sizeoflen);
      len[0]=lentot;
    }
    /* if array is not small, then all elements are distributed evenly across the compute processes */
    else {
      int minlen = lentot / size;
      int nbiglen= lentot % size;

      if (minlen==0) sizeoflen=lentot;
      else sizeoflen=size;

      len.resize(sizeoflen);

      for(int i=0; i<sizeoflen; i++) {
         if(i<nbiglen) len[i]=minlen+1;
         else len[i]=minlen;
      }
    }

    twosided_helpga_create_irreg(mproc, len.data(), sizeoflen, handle, name, dtype);

    if (DEBUG_) printf("%5d: twosided_helpga_create: end. mproc=%d, handle=%d\n",ProcID(),mproc,*handle);

    return 0;

}


int twosided_helpga_locate_server( int handle, int ilo, int ihigh, int *map, int *proclist_help, int *np_help)
{
    int i, rank;
    int size;
    int lenleft,iilow,iihig,offset;
    int iproclow=0,iprochigh=0;
    int handle_orig;
    int lentot,*len_help=NULL;
    MPIHELPGA helpga=NULL;

    if (DEBUG_) printf("%5d: In twosided_helpga_locate_server: begin. handle=%d, ilo=%d, ihigh=%d\n",ProcID(),handle,ilo,ihigh);

    handle_orig=twosided_helpga_handle_orig(handle);
    helpga=twosided_helpga_index[handle_orig].ptr;

    size=NUMBER_OF_SERVER;
    lentot=helpga->nele;
    len_help=helpga->len_help;

    if ( ilo < 1 || ilo >ihigh ||ihigh > lentot) {
       fprintf(stderr,"ERROR in twosided_helpga_locate_server: over range! ilo=%d,ihigh=%d\n",ilo,ihigh);
       MPI_Abort(MPI_COMM_WORLD,MPI_ERR_OTHER);
    }

    for (lenleft=0,i=0;i<size;i++){
       if (ilo   >=lenleft+1 && ilo   <= lenleft+len_help[i]) iproclow=i;
       if (ihigh >=lenleft+1 && ihigh <= lenleft+len_help[i]) { iprochigh=i; break; }
       lenleft=lenleft + len_help[i];
    }

    *np_help=iprochigh-iproclow+1;

    for (lenleft=0,i=0;i<iproclow;i++) lenleft=lenleft + len_help[i];

    for (rank=iproclow;rank<=iprochigh;rank++) {
       i=rank-iproclow;
       proclist_help[i]=RankNumber_of_Server(rank);
       iilow=lenleft+1;
       iihig=lenleft + len_help[rank];
       offset=2*i;
       map[offset]  = iilow < ilo   ? ilo:   iilow;
       map[offset+1]= iihig > ihigh ? ihigh: iihig;
       lenleft=iihig;
       if (DEBUG_) printf("%5d: In twosided_helpga_locate_server: i=%d, proclist_help[i]=%d, map=%d  %d\n",ProcID(),i,proclist_help[i],map[offset],map[offset+1]);
    }

    len_help=NULL;
    helpga=NULL;
    if (DEBUG_) printf("%5d: In twosided_helpga_locate_server: end. handle=%d, ilo=%d, ihigh=%d\n",ProcID(),handle,ilo,ihigh);

    return 0;
}

/* find out the range of helpga that process iproc "owns". iproc should be serial rank of any valid compute process */
/* If no array elements are owned by process iproc, the range is returned as ilo[ ]=0 and ihigh[ ]= -1  */
/* Note: here the process rank is the serial rank of compute process */
int twosided_helpga_distrib( int handle, int iproc, int *ilo, int *ihigh)
{
    int i;
    int size;
    int lenleft;
    int handle_orig;
    MPIHELPGA helpga=NULL;

    if (DEBUG_) printf("%5d: In twosided_helpga_distrib begin: handle=%d\n",ProcID(),handle);

    handle_orig=twosided_helpga_handle_orig(handle);
    helpga=twosided_helpga_index[handle_orig].ptr;

    size=NProcs_Work();

    if ( iproc < 0 || iproc >=size ) {
       fprintf(stderr,"ERROR in twosided_helpga_distrib: iproc= %d over range!!\n",iproc);
       MPI_Abort(MPI_COMM_WORLD,MPI_ERR_OTHER);
    }

    for (lenleft=0,i=0;i<iproc;i++) lenleft=lenleft + helpga->len[i];

    if (helpga->len[iproc] !=0 ) {
       *ilo=lenleft+1;
       *ihigh=lenleft+helpga->len[iproc];
    }
    else {
       *ilo=0;
       *ihigh=-1;
    }

    helpga=NULL;

    if (DEBUG_) printf("%5d: In twosided_helpga_distrib end: handle=%d, iproc=%d, ilo=%d, ihigh=%d\n",ProcID(),handle,iproc,*ilo,*ihigh);

    return 0;
}


/* Return the list of the helpga compute processes id that 'own' the data. Parts of the specified patch might be actually 'owned' by several processes.
 If lo/hi are out of bounds then error is given. np is equal to the number of processes that hold the data . */
/* Note: here the process rank is the serial rank of compute process */
int twosided_helpga_location( int handle, int ilo, int ihigh, int *map, int *proclist, int *np)
{
    int i, rank;
    int size;
    int lenleft,iilow,iihig,offset;
    int iproclow=0,iprochigh=0;
    int lentot;
    int handle_orig;
    MPIHELPGA helpga=NULL;

    if (DEBUG_) printf("%5d: In twosided_helpga_location: begin. handle=%d\n",ProcID(),handle);

    handle_orig=twosided_helpga_handle_orig(handle);
    helpga=twosided_helpga_index[handle_orig].ptr;

    size=NProcs_Work();
    lentot=helpga->nele;

    if ( ilo < 1 || ilo >ihigh ||ihigh > lentot) {
       fprintf(stderr,"ERROR in twosided_helpga_location: over range! ilo=%d,ihigh=%d\n",ilo,ihigh);
       MPI_Abort(MPI_COMM_WORLD,MPI_ERR_OTHER);
    }
    for (lenleft=0,i=0;i<size;i++){
       if (ilo   >=lenleft+1 && ilo   <= lenleft+helpga->len[i]) iproclow=i;
       if (ihigh >=lenleft+1 && ihigh <= lenleft+helpga->len[i]) { iprochigh=i; break; }
       lenleft=lenleft + helpga->len[i];
    }

    *np=iprochigh-iproclow+1;

    for (lenleft=0,i=0;i<iproclow;i++) lenleft=lenleft + helpga->len[i];

    for (rank=iproclow;rank<=iprochigh;rank++) {
       i=rank-iproclow;
       proclist[i]=rank;
       iilow=lenleft+1;
       iihig=lenleft + helpga->len[rank];
       offset=2*i;
       map[offset]  = iilow < ilo   ? ilo:   iilow;
       map[offset+1]= iihig > ihigh ? ihigh: iihig;
       lenleft=iihig;
    }

    helpga=NULL;

    if (DEBUG_) printf("%5d: In twosided_helpga_location end. handle=%d, ilo=%d, ihigh=%d, np=%d\n",ProcID(),handle,ilo,ihigh,*np);

    return 0;
}

/* one-element rma operations for helpga */
/* MPI_Datatye can only be integer type (MPI_INT32_T and MPI_INT64_T but not MPI_DOUBLE) */
int64_t twosided_helpga_one(int mproc, int64_t nelem_valput, int ielem, int handle)
/*
  Operations for helpga:
________________________________________________________________________________________________________
|                |     COLLECFLAG       |       RMAONEFLAG             |        RMAETRFLAG             |
|  mproc = 0 ... |   create helpga      | get one element from helpga  | get elements from helpga      |
|  mproc > 0 ... |   zeroize helpga     | fetch-and-add to helpga      | accumulate helpga elements    |
|  mproc < 0 ... |   destroy helpga     | put one value to helpga      | put values to helpga elements |
--------------------------------------------------------------------------------------------------------
*/
{
  int  type = RMAONEFLAG;
  dataserver buf[4];
  int64_t local=0;
  MPI_Status status;
  int  handle_orig;
  MPI_Datatype dtype;             /* MPI Datatype */
  int  server=0;                  /* id of server process */
  int  lentot;                    /* total number of elements     */
  int  size,lenleft,i,iserver=0;
  int  *len_help=NULL;
  MPIHELPGA helpga=NULL;

  if (DEBUG_) printf("%5d: twosided_helpga_one: begin. type=%d, mproc=%d, handle=%d, ielem=%d\n",ProcID(),type,mproc,handle,ielem);

  handle_orig=twosided_helpga_handle_orig(handle);
  helpga=twosided_helpga_index[handle_orig].ptr;

  dtype=helpga->dtype;
  lentot=helpga->nele;
  len_help=helpga->len_help;

  if (dtype==MPI_INT32_T || dtype==MPI_INT64_T) {
    if (DEBUG_) printf("%5d: twosided_helpga_one: array with handle=%d  is an integer MPIGA.\n",ProcID(),handle);
  }
  else  MPIGA_Error(" twosided_helpga_one: wrong MPI_Datatype ",0);
  if ( ielem < 1 || ielem > lentot) {
    fprintf(stderr,"  twosided_helpga_one ERROR: over range! lentot=%d, ielem=%d\n",lentot,ielem);
    MPIGA_Error(" twosided_helpga_one: ielem over range ",0);
  }

  if (DEBUG_) {
    printf("%5d: twosided_helpga_one: mid. mproc=%d, handle=%d\n",ProcID(),mproc,handle);
    fflush(stdout);
  }

  if (SR_parallel) {
     buf[0] = (dataserver)mproc; /* COLLECFLAG: create(=0), zeroize(>0), destroy(<0); RMAONEFLAG: get(=0), fetch-and-add(>0), put(<0) */
     buf[1] = (dataserver)nelem_valput; /* COLLECFLAG and mproc=0: number of elements; RMAONEFLAG: value to be put(mproc<0), increment value(mproc>0); others: no use */
     buf[2] = (dataserver)ielem;        /* COLLECFLAG: MPI_Datatype; RMAONEFLAG: sequence number of element (1,2,...,n) */
     buf[3] = (dataserver)handle;      /* sequence number of helpga */

/* find out in which helper process the element is located */
     size=NUMBER_OF_SERVER;
     for(lenleft=0,i=0;i<size;i++){
       if (ielem >= lenleft+1 && ielem  <= lenleft+len_help[i]) { iserver=i; break; }
       lenleft=lenleft + len_help[i];
     }
     server=RankNumber_of_Server(iserver);

     buf[2] = (dataserver)(ielem-lenleft); /* RMAONEFLAG: sequence number of element (1,2,...,len_help[i]) in len_help */

     if ( buf[2] >  (dataserver) len_help[iserver] || buf[2] <= (dataserver)0 ) {
       printf("%5d: twosided_helpga_one: ERROR mproc=%d, handle=%d,ielem_inhelp=%ld,len_help[i]=%d\n",ProcID(),mproc,handle,(long)buf[2],len_help[iserver]);
       MPIGA_Error(" twosided_helpga_one: overange ",0);
     }

     if (use_helper_server) {
       MPI_Send(buf, 4, DATASERVER_MPI,  server, type, MPI_COMM_WORLD);
       MPI_Recv(buf, 1, DATASERVER_MPI,  server, type, MPI_COMM_WORLD, &status);
       local=(int64_t)buf[0];
     }
  }
  else {
     /* Not running in parallel ... just do a simulation */
    int32_t *i32buf=NULL;
    int64_t *i64buf=NULL;
    if (mproc == 0) {                 /* get a value from local helpga      */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)helpga->ptr_buf;
         local = (int64_t)i32buf[ielem-1];
         i32buf=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)helpga->ptr_buf;
         local = (int64_t)i64buf[ielem-1];
         i64buf=NULL;
      }
    }
    else if (mproc == 1) {            /* fetch-and-add INCR to local helpga */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)helpga->ptr_buf;
         local = (int64_t)i32buf[ielem-1];
         i32buf[ielem-1]+=(int32_t)nelem_valput;
         i32buf=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)helpga->ptr_buf;
         local = (int64_t)i64buf[ielem-1];
         i64buf[ielem-1]+=(int64_t)nelem_valput;
         i64buf=NULL;
      }
    }
    else if (mproc == -1) {           /* put a value  to local helpga       */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)helpga->ptr_buf;
         local = (int64_t)i32buf[ielem-1];
         i32buf[ielem-1]=(int32_t)nelem_valput;
         i32buf=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)helpga->ptr_buf;
         local = (int64_t)i64buf[ielem-1];
         i64buf[ielem-1]=(int64_t)nelem_valput;
         i64buf=NULL;
      }
    }
    else  {
      MPIGA_Error("twosided_helpga_one: sequential version with silly mproc ", mproc);
    }
  } /* sequential case */
  len_help=NULL;
  helpga=NULL;

  if (DEBUG_) printf("%5d: twosided_helpga_one: end. type=%d, mproc=%d, handle=%d, ielem=%d\n",ProcID(),type,mproc,handle,ielem);

  return local;
}


/* many-element rma operations for helpga, MPI_Datatye can be MPI_INT32_T, MPI_INT32_T or MPI_DOUBLE */
void twosided_helpga_extra(int mproc, int nelem, int ielem, int handle, void *buff)
/*
  Operations for helpga:
________________________________________________________________________________________________________
|                |     COLLECFLAG       |       RMAONEFLAG             |        RMAETRFLAG             |
|  mproc = 0 ... |   create helpga      | get one element from helpga  | get elements from helpga      |
|  mproc > 0 ... |   zeroize helpga     | fetch-and-add to helpga      | accumulate helpga elements    |
|  mproc < 0 ... |   destroy helpga     | put one value to helpga      | put values to helpga elements |
--------------------------------------------------------------------------------------------------------
*/
{
  int  type = RMAETRFLAG;
  dataserver buf[4];
  int  handle_orig;
  MPI_Datatype dtype;            /* MPI Datatype */
  int sizeofdtype;               /* Size of MPI Datatype */
  int  server;                   /* id of server process */
  int  np_help, lenleft,lenleft_save;
  int  i,ilo,ihigh,ifirst,ilast,iserver_first,iserver;
  int  *len_help=NULL;
  MPIHELPGA helpga=NULL;

  if (DEBUG_) printf("%5d: twosided_helpga_extra: begin. type=%d, mproc=%d, handle=%d\n",ProcID(),type,mproc,handle);

  handle_orig=twosided_helpga_handle_orig(handle);
  helpga=twosided_helpga_index[handle_orig].ptr;

  dtype=helpga->dtype;
  len_help=helpga->len_help;

  MPI_Type_size( dtype, &sizeofdtype );

  if (SR_parallel) {
   if (use_helper_server) {

     buf[0] = (dataserver)mproc; /* COLLECFLAG: create(=0), zeroize(>0), destroy(<0); RMAONEFLAG: get(=0), fetch-and-add(>0), put(<0); RMAETRFLAG: get n elements(=0), put n elements(<0) */
     buf[1] = (dataserver)nelem; /* RMAETRFLAG: number of elements to be gotten/put/accumulated */
     buf[2] = (dataserver)ielem; /* RMAETRFLAG: sequence number of element (1,2,...,n) */
     buf[3] = (dataserver)handle;     /* sequence number of helpga */

     ilo=ielem;
     ihigh=ielem+nelem-1;
     twosided_helpga_locate_server(handle, ilo, ihigh, twosided_helpga_map, twosided_helpga_proclist, &np_help);
/* iserver_first: serial number of server for lowest element [ilo]; lenleft: number of elements stored in helper servers < iserver_first */
     iserver_first=SerialNumber_of_Server(twosided_helpga_proclist[0]);
     for (lenleft=0,i=0;i<iserver_first;i++) lenleft=lenleft + len_help[i];

     if (mproc > 0 || mproc < 0) {
/* mproc>0, accumulate a set of elements to helpga, here ensure it is an atomic operation */
/* mproc<0, put a set of elements to helpga */
      lenleft_save=lenleft;
      for(i=0;i<np_help;i++) {
       server = twosided_helpga_proclist[i];
       ifirst = twosided_helpga_map[2*i];
       ilast  = twosided_helpga_map[2*i+1];

       ielem = ifirst-lenleft;
       nelem = ilast-ifirst+1;

       buf[1] = (dataserver)nelem; /* RMAETRFLAG: number of elements to be gotten/put */
       buf[2] = (dataserver)ielem; /* RMAETRFLAG: sequence number of element (1,2,...,n) */

       MPI_Send(buf, 4, DATASERVER_MPI,  server, type, MPI_COMM_WORLD);
       iserver=SerialNumber_of_Server(server);
       lenleft=lenleft + len_help[iserver];
      }
      lenleft=lenleft_save;
     }

     std::vector<MPI_Request> requests2(np_help);
     std::vector<MPI_Request> requests3(np_help);

/* get/put/accumulate the data from/to servers */
     for(i=0;i<np_help;i++) {
       server = twosided_helpga_proclist[i];
       ifirst = twosided_helpga_map[2*i];
       ilast  = twosided_helpga_map[2*i+1];

       ielem = ifirst-lenleft;
       nelem = ilast-ifirst+1;

       buf[1] = (dataserver)nelem; /* RMAETRFLAG: number of elements to be gotten/put */
       buf[2] = (dataserver)ielem; /* RMAETRFLAG: sequence number of element (1,2,...,n) */

       if (mproc == 0) {                 /* get a set of elements from helpga  */
         MPI_Send(buf, 4, DATASERVER_MPI,  server, type, MPI_COMM_WORLD);

         MPI_Irecv(buff, nelem, dtype, server, type, MPI_COMM_WORLD, &requests2[i]);
         if (DEBUG_) {
           MPI_Recv(NULL, 0, MPI_BYTE, server, WAKEUPTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* SUCCESS: receive for notification from server */
         }
       }
       else  if (mproc < 0) {            /* put a set of elements to helpga */
         MPI_Isend(buff, nelem, dtype, server, type, MPI_COMM_WORLD, &requests2[i]);
         MPI_Irecv(NULL, 0, MPI_BYTE, server, WAKEUPTAG, MPI_COMM_WORLD, &requests3[i]); /* SUCCESS: receive for notification from server */
       }
       else { /* mproc>0, accumulate a set of elements to helpga */
         MPI_Isend(buff, nelem, dtype, server, type, MPI_COMM_WORLD, &requests2[i]);
         MPI_Irecv(NULL, 0, MPI_BYTE, server, WAKEUPTAG, MPI_COMM_WORLD, &requests3[i]); /* SUCCESS: receive for notification from server */
       }

       iserver=SerialNumber_of_Server(server);
       lenleft=lenleft + len_help[iserver];
       buff = (void *)( ((char *)buff) + nelem * sizeofdtype );
     }

     if (mproc == 0) { /* get */
       MPI_Waitall(np_help, requests2.data(), MPI_STATUSES_IGNORE);
     }
     else { /* put and accumulate */
       MPI_Waitall(np_help, requests2.data(), MPI_STATUSES_IGNORE);
       MPI_Waitall(np_help, requests3.data(), MPI_STATUSES_IGNORE);
     }
   } /* end of use_helper_server loop */
  }
  else {
     /* Not running in parallel ... just do a simulation */
    int32_t *i32buf=NULL,*i32buf_helpga=NULL;
    int64_t *i64buf=NULL,*i64buf_helpga=NULL;
    double *dbuf=NULL,*dbuf_helpga=NULL;
    if (mproc == 0) {                 /* get a set of values from local helpga      */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)buff;
         i32buf_helpga=(int32_t *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) i32buf[i]=i32buf_helpga[ielem-1+i];
         i32buf=NULL;
         i32buf_helpga=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)buff;
         i64buf_helpga=(int64_t *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) i64buf[i]=i64buf_helpga[ielem-1+i];
         i64buf=NULL;
         i64buf_helpga=NULL;
      }
      else if (dtype==MPI_DOUBLE) {
         dbuf=(double *)buff;
         dbuf_helpga=(double *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) dbuf[i]=dbuf_helpga[ielem-1+i];
         dbuf=NULL;
         dbuf_helpga=NULL;
      }
      else {
         MPIGA_Error("twosided_helpga_extra [get]: wrong MPI_Datatype ",0);
      }
    }
    else if (mproc == -1) {           /* put a set of values  to local helpga       */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)buff;
         i32buf_helpga=(int32_t *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) i32buf_helpga[ielem-1+i]=i32buf[i];
         i32buf=NULL;
         i32buf_helpga=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)buff;
         i64buf_helpga=(int64_t *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) i64buf_helpga[ielem-1+i]=i64buf[i];
         i64buf=NULL;
         i64buf_helpga=NULL;
      }
      else if (dtype==MPI_DOUBLE) {
         dbuf=(double *)buff;
         dbuf_helpga=(double *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) dbuf_helpga[ielem-1+i]=dbuf[i];
         dbuf=NULL;
         dbuf_helpga=NULL;
      }
      else {
         MPIGA_Error("twosided_helpga_extra [put]: wrong MPI_Datatype ",0);
      }
    }
    else if (mproc == 1) {  /* mproc>0, accumulate a set of elements to local helpga */
      if (dtype==MPI_INT32_T) {
         i32buf=(int32_t *)buff;
         i32buf_helpga=(int32_t *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) i32buf_helpga[ielem-1+i]+=i32buf[i];
         i32buf=NULL;
         i32buf_helpga=NULL;
      }
      else if (dtype==MPI_INT64_T) {
         i64buf=(int64_t *)buff;
         i64buf_helpga=(int64_t *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) i64buf_helpga[ielem-1+i]+=i64buf[i];
         i64buf=NULL;
         i64buf_helpga=NULL;
      }
      else if (dtype==MPI_DOUBLE) {
         dbuf=(double *)buff;
         dbuf_helpga=(double *)helpga->ptr_buf;
         for (i=0;i<nelem;i++) dbuf_helpga[ielem-1+i]+=dbuf[i];
         dbuf=NULL;
         dbuf_helpga=NULL;
      }
      else {
         MPIGA_Error("twosided_helpga_extra [acc]: wrong MPI_Datatype ",0);
      }
    }
    else {
      MPIGA_Error("twosided_helpga_extra: sequential version with silly mproc ", mproc);
    }
  } /* sequential case */
  len_help=NULL;
  helpga=NULL;

  if (DEBUG_) printf("%5d: twosided_helpga_extra: end. type=%d, mproc=%d, handle=%d\n",ProcID(),type,mproc,handle);
}

/* many-element accumulation operation for helpga, MPI_Datatye can be MPI_INT32_T, MPI_INT32_T or MPI_DOUBLE */
/* determine if fac==1/1.0d0 */
void twosided_helpga_extra_acc(int mproc, int nelem, int ielem, int handle, void *buf, void *fac)
{
    MPI_Datatype dtype;  /* MPI Datatype for helpga element */
    void *alphabuf=NULL;
    int32_t *i32tempbuf=NULL,*i32fac=NULL;
    int64_t *i64tempbuf=NULL,*i64fac=NULL;
    double *dtempbuf=NULL,*dfac=NULL;
    std::vector<int32_t> i32alphabuf;
    std::vector<int64_t> i64alphabuf;
    std::vector<double> dalphabuf;

    if (DEBUG_) printf("%5d: twosided_helpga_extra_acc: begin. handle=%d\n",ProcID(),handle);
    dtype=twosided_helpga_inquire_dtype(handle);
    if (dtype==MPI_INT32_T) {
       i32fac=(int32_t *)fac;
       if ((*i32fac)==(int32_t)1) alphabuf=buf;
       else {
          i32tempbuf=(int32_t *)buf;
	  i32alphabuf.resize(nelem);
          for (int i=0;i<nelem;i++)i32alphabuf[i]=(*i32fac)*i32tempbuf[i];
          alphabuf=(void *)i32alphabuf.data();
       }
    }
    else if (dtype==MPI_INT64_T) {
       i64fac=(int64_t *)fac;
       if ((*i64fac)==(int64_t)1) alphabuf=buf;
       else {
          i64tempbuf=(int64_t *)buf;
	  i64alphabuf.resize(nelem);
          for (int i=0;i<nelem;i++)i64alphabuf[i]=(*i64fac)*i64tempbuf[i];
          alphabuf=(void *)i64alphabuf.data();
       }
    }
    else if (dtype==MPI_DOUBLE) {
       dfac=(double *)fac;
       if (std::abs((*dfac)-1.0e0)<1.0e-6) alphabuf=buf;
       else {
          dtempbuf=(double *)buf;
	  dalphabuf.resize(nelem);
          for (int i=0;i<nelem;i++) dalphabuf[i]=(*dfac)*dtempbuf[i];
          alphabuf=(void *)dalphabuf.data();
       }
    }
    else {
       MPIGA_Error("twosided_helpga_extra_acc: wrong MPI_Datatype ",0);
    }

    twosided_helpga_extra(mproc, nelem, ielem, handle, alphabuf);

    if (DEBUG_) printf("%5d: twosided_helpga_extra_acc: end. handle=%d\n",ProcID(),handle);
}


/* get the original handle of helpga, and check whether handle is out of range or allocated */
int twosided_helpga_handle_orig( int handle )
{
    int handle_orig;

    if (DEBUG_) printf("%5d: twosided_helpga_handle_orig: begin. handle=%d\n",ProcID(),handle);

    /* check whether helpga handle is out of range, and check whether it is allocated */
    handle_orig=handle-TWOSIDED_HELPGA_OFFSET;
    if(handle_orig < 0 || handle_orig >= MAX_TWOSIDED_HELPGA_ARRAYS){
       fprintf(stderr," twosided_helpga_handle_orig ERROR:  invalid handle [%d]. Should be [ %d -- %d ].\n",
                       handle,TWOSIDED_HELPGA_OFFSET,(TWOSIDED_HELPGA_OFFSET+MAX_TWOSIDED_HELPGA_ARRAYS));
       exit (1);
    }
    if(twosided_helpga_index[handle_orig].actv==0){
       fprintf(stderr," twosided_helpga_handle_orig ERROR:  no memory allocated for this helpga [%d].\n",handle);
       exit (1);
    }
    if (DEBUG_) printf("%5d: twosided_helpga_handle_orig: end. handle=%d\n",ProcID(),handle);
    return handle_orig;
}


/* get the MPI_Datatype of a helpga represented by handle */
MPI_Datatype twosided_helpga_inquire_dtype( int handle )
{
    int handle_orig;
    MPI_Datatype mpidtype;
    MPIHELPGA helpga=NULL;

    if (DEBUG_) printf("%5d: twosided_helpga_inquire_dtype: begin. handle=%d\n",ProcID(),handle);

    handle_orig=twosided_helpga_handle_orig(handle);
    helpga=twosided_helpga_index[handle_orig].ptr;

    mpidtype=helpga->dtype;

    helpga=NULL;

    if (DEBUG_) printf("%5d: twosided_helpga_inquire_dtype: end. handle=%d\n",ProcID(),handle);

    return mpidtype;
}


/* get the name of a helpga represented by handle */
int twosided_helpga_inquire_name( int handle, char **name )
{
    int handle_orig;
    MPIHELPGA helpga=NULL;

    if (DEBUG_) printf("%5d: twosided_helpga_inquire_name: begin. handle=%d\n",ProcID(),handle);

    handle_orig=twosided_helpga_handle_orig(handle);
    helpga=twosided_helpga_index[handle_orig].ptr;

    *name=helpga->name;

    helpga=NULL;

    if (DEBUG_) printf("%5d: twosided_helpga_inquire_name: end. handle=%d, name=%s\n",ProcID(),handle,*name);
    return 0;
}


/* Collective operations(allocate and free) for mutexes */
/* Creates and destroys a set of mutexes (in total: number). Only one set of mutexes can exist at a time. Mutexes can be
created and destroyed as many times as needed. Mutexes are numbered: 0, ..., number-1. */
int twosided_helpmutex_collect(int mproc, int inum, int number)
/*
  Operations for helpmutex:
__________________________________________________________________
|                |     MUTCOLFLAG       |      MUTLOCFLAG        |
|  mproc > 0 ... | create helpmutexes   |   lock a helpmutex     |
|  mproc < 0 ... | destroy helpmutexes  | unlock a helpmutex     |
------------------------------------------------------------------
MUTCOLFLAG(inum): first sequential number of mutexes to be allocated (mproc > 0) or freed (mproc < 0)
MUTCOLFLAG(number): number of mutexes to be allocated (mproc > 0) or freed (mproc < 0)
MUTLOCFLAG(inum): sequential number of mutex to be locked (mproc > 0) and unlocked (mproc < 0)
*/
{
  int  type = MUTCOLFLAG;
  dataserver buf[3];
  int  local=0;
  int  server;                    /* id of server process */
  int  i;

  if (DEBUG_) printf("%5d: twosided_helpmutex_collect: begin. type=%d, mproc=%d, inum=%d, number=%d\n",ProcID(),type,mproc,inum,number);

  if ( MPIGA_WORK_COMM != MPI_COMM_NULL) MPI_Barrier(MPIGA_WORK_COMM);

  if ( number <=0 || number > MAX_TWOSIDED_HELPMUTEX_ARRAYS ) {
    fprintf(stderr," twosided_helpmutex_lock ERROR:  over range! number= %d,should be [ 1 -- %d ].\n",number,MAX_TWOSIDED_HELPMUTEX_ARRAYS);
    exit (1);
  }

  if( mproc==0 ) { /* should be >0(create) or <0(destroy), can't =0 */
    fprintf(stderr," twosided_helpmutex_collect ERROR:  not create/destroy operation.\n");
    exit (1);
  }
  else if ( mproc>0 ) { /* allocate operation */
    for(i=inum;i<inum+number; i++) {
       twosided_helpmutex_index[i].actv = 1;
    }
    if (inum>=MAX_TWOSIDED_HELPMUTEX_GA) twosided_helpmutex_num=number;
  }
  else { /* free operation */
    for(i=inum;i<inum+number; i++) {
       twosided_helpmutex_index[i].actv = 0;
       if ( twosided_helpmutex_index[i].mutex != NULL ) { free (twosided_helpmutex_index[i].mutex); twosided_helpmutex_index[i].mutex=NULL; }
    }
    if (inum>=MAX_TWOSIDED_HELPMUTEX_GA) twosided_helpmutex_num=0;
  }

  if (SR_parallel) {
     buf[0] = (dataserver)mproc;   /* MUTCOLFLAG: allocate(mproc>0), free(mproc<0);  MUTLOCFLAG: lock(mproc>0), unlock(mproc<0). */
     buf[1] = (dataserver)inum;    /* MUTCOLFLAG(inum): first sequential number of mutexes to be allocated (mproc > 0) or freed (mproc < 0)  */
     buf[2] = (dataserver)number;  /* MUTCOLFLAG(number): number of mutexes to be allocated (mproc > 0) or freed (mproc < 0) */

     if (use_helper_server) {
       server = LastServerID();                     /* helpmutex server is always the last process(ie, NXTVAL server) */
       MPI_Send(buf, 3, DATASERVER_MPI, server, type, MPI_COMM_WORLD);
       MPI_Recv(NULL, 0, MPI_BYTE, server, type, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* SUCCESS: receive notification from server */
     }
  }
  else {
     /* sequential case: Not running in parallel ... just do a simulation */
     if ( abs(mproc) != 1 )  MPIGA_Error("twosided_helpmutex_collect: sequential version with silly mproc ", mproc);
     local = 0;
  } /* sequential case */


  if ( MPIGA_WORK_COMM != MPI_COMM_NULL) MPI_Barrier(MPIGA_WORK_COMM);

  if (DEBUG_) printf("%5d: twosided_helpmutex_collect: end. type=%d, mproc=%d, inum=%d, number=%d\n",ProcID(),type,mproc,inum,number);

  return local;
}


/* Lock and unlock operations for a specific mutex */
/* Lock/unlock a mutex object identified by the original mutex number. It is a fatal error for a process
   to attempt to lock/unlock a mutex which has already been locked/unlocked by this process */
int twosided_helpmutex_lock(int mproc, int inum)
/*
  Operations for helpmutex:
__________________________________________________________________
|                |     MUTCOLFLAG       |      MUTLOCFLAG        |
|  mproc > 0 ... | create helpmutexes   |   lock a helpmutex     |
|  mproc < 0 ... | destroy helpmutexes  | unlock a helpmutex     |
------------------------------------------------------------------
MUTCOLFLAG(inum): first sequential number of mutexes to be allocated (mproc > 0) or freed (mproc < 0)
MUTCOLFLAG(number): number of mutexes to be allocated (mproc > 0) or freed (mproc < 0)
MUTLOCFLAG(inum): sequential number of mutex to be locked (mproc > 0) and unlocked (mproc < 0)
*/
{
  int  type = MUTLOCFLAG;
  dataserver buf[2];
  int  local=0;
  int  server;                    /* id of server process */

  if (DEBUG_) printf("%5d: twosided_helpmutex_lock: begin. type=%d, mproc=%d, inum(mutex)=%d\n",ProcID(),type,mproc,inum);

  if ( inum <0 || inum >= MAX_TWOSIDED_HELPMUTEX_ARRAYS ) {
    fprintf(stderr," twosided_helpmutex_lock ERROR:  over range! original inum(mutex)=%d, max num=%d\n",inum,MAX_TWOSIDED_HELPMUTEX_ARRAYS);
    exit (1);
  }
  if ( twosided_helpmutex_index[inum].actv==0 ) {
    fprintf(stderr," twosided_helpmutex_lock ERROR: original inum(mutex)=%d has not been activated yet.\n",inum);
    exit (1);
  }

  if( mproc==0 ) { /* should be >0(lock) or <0(unlock), can't =0 */
    fprintf(stderr," twosided_helpmutex_lock ERROR:  not lock/unlock operation.\n");
    exit (1);
  }
  else if ( mproc>0 ) { /* lock operation */
    if ( twosided_helpmutex_index[inum].lock==1 ) {
       fprintf(stderr," twosided_helpmutex_lock ERROR:  attempt to lock a mutex [%d] which has already been locked.\n",inum);
       exit (1);
    }
    else {
       twosided_helpmutex_index[inum].lock = 1;
    }
  }
  else { /* unlock operation */
    if ( twosided_helpmutex_index[inum].lock==0 ) {
       fprintf(stderr," twosided_helpmutex_lock ERROR:  attempt to unlock a mutex [%d] which has not been locked.\n",inum);
       exit (1);
    }
    else {
       twosided_helpmutex_index[inum].lock = 0;
    }
  }

  if (SR_parallel) {
     buf[0] = (dataserver)mproc;   /* MUTCOLFLAG: allocate(mproc>0), free(mproc<0);  MUTLOCFLAG: lock(mproc>0), unlock(mproc<0). */
     buf[1] = (dataserver)inum;    /* MUTLOCFLAG(inum): sequential number of mutex to be locked (mproc>0) and unlocked (mproc<0) */

     if (use_helper_server) {
       server = LastServerID();                     /* helpmutex server is always the last process(ie, NXTVAL server) */
       MPI_Send(buf, 2, DATASERVER_MPI, server, type, MPI_COMM_WORLD);
       MPI_Recv(NULL, 0, MPI_BYTE, server, type, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* SUCCESS: receive notification from server */
     }
  }
  else {
     /* sequential case: Not running in parallel ... just do a simulation */
     if ( abs(mproc) != 1 )  MPIGA_Error("twosided_helpmutex_lock: sequential version with silly mproc ", mproc);
     local = 0;
  } /* sequential case */


  if (DEBUG_) printf("%5d: twosided_helpmutex_lock: end. type=%d, mproc=%d, inum(mutex)=%d\n",ProcID(),type,mproc,inum);

  return local;
}

 /* lock a twosided_helpmutex object identified by the original mutex number. */
int lock_twosided_helpmutex_orig(int inum_orig)
{
    int mproc=NProcs_Work();
    int mpierr=0;

    mpierr=twosided_helpmutex_lock(mproc, inum_orig);
    return mpierr;
}

 /* unlock a twosided_helpmutex object identified by the original mutex number. */
int unlock_twosided_helpmutex_orig(int inum_orig)
{
    int mproc=-NProcs_Work();
    int mpierr;

    mpierr=twosided_helpmutex_lock(mproc, inum_orig);
    return mpierr;
}

/* lock a mutex object identified by the wrapped mutex number. It is a fatal error for a process
   to attempt to lock a mutex which has already been locked by this process */
/* It is only called by NON-GA routines, CAN NOT be called inside mpiga operations */
int lock_twosided_helpmutex(int inum)
{
    int mproc=NProcs_Work();
    int inum_orig;

    if (DEBUG_) printf("%5d: In lock_twosided_helpmutex begin: mutex num=%d, total num=%d\n",ProcID(),inum,twosided_helpmutex_num);

    if (inum <0 || inum >= twosided_helpmutex_num ) {
       fprintf(stderr," lock_twosided_helpmutex ERROR: over range! mutex num=%d, total num=%d\n",inum,twosided_helpmutex_num);
       exit (1);
    }
    else {
       inum_orig=inum+MAX_TWOSIDED_HELPMUTEX_GA;
       twosided_helpmutex_lock(mproc, inum_orig);
    }

    if (DEBUG_) printf("%5d: In lock_twosided_helpmutex end: mutex %d has been locked.\n",ProcID(),inum);

    return 0;
}

/* unlock a mutex object identified by the wrapped mutex number. It is a fatal error for a process
 * to attempt to unlock a mutex which has not been locked by this process. */
/* It is only called by NON-GA routines, CAN NOT be called inside mpiga operations */
int unlock_twosided_helpmutex(int inum)
{
    int mproc=-NProcs_Work();
    int inum_orig;

    if (DEBUG_) printf("%5d: In lock_twosided_helpmutex begin: mutex num=%d, total num=%d\n",ProcID(),inum,twosided_helpmutex_num);

    if (inum <0 || inum >= twosided_helpmutex_num ) {
       fprintf(stderr," lock_twosided_helpmutex ERROR: over range! mutex num=%d, total num=%d\n",inum,twosided_helpmutex_num);
       exit (1);
    }
    else {
       inum_orig=inum+MAX_TWOSIDED_HELPMUTEX_GA;
       twosided_helpmutex_lock(mproc, inum_orig);
    }

    if (DEBUG_) printf("%5d: In lock_twosided_helpmutex end: mutex %d has been unlocked.\n",ProcID(),inum);

    return 0;
}

/* allocate one mutex object identified by the original mutex number. Check whether it has overranged, activated, or locked */
int alloc_twosided_helpmutex_orig(int inum_orig)
{
    int mproc=NProcs_Work();
    int number=1;

    if (DEBUG_) printf("%5d: In alloc_twosided_helpmutex_orig begin: original mutex num=%d\n",ProcID(),inum_orig);

    if (inum_orig <0 || inum_orig >= MAX_TWOSIDED_HELPMUTEX_ARRAYS ) {
       fprintf(stderr,"ERROR in alloc_twosided_helpmutex_orig: over range! original mutex num=%d, max num=%d\n",inum_orig,MAX_TWOSIDED_HELPMUTEX_ARRAYS);
       exit (1);
    }
    else if ( twosided_helpmutex_index[inum_orig].actv==1 ) {
       fprintf(stderr,"ERROR in alloc_twosided_helpmutex_orig: original mutex %d has already been activated.\n",inum_orig);
       exit (1);
    }
    else if ( twosided_helpmutex_index[inum_orig].lock==1 ) {
       fprintf(stderr,"ERROR in alloc_twosided_helpmutex_orig: original mutex %d has been locked before allocating.\n",inum_orig);
       exit (1);
    }
    else {
       twosided_helpmutex_collect(mproc, inum_orig, number);
    }

    if (DEBUG_) printf("%5d: In alloc_twosided_helpmutex_orig end: original mutex=%d\n",ProcID(),inum_orig);

    return 0;
}

/* free one mutex object identified by the original mutex number. Check whether it has overranged, activated, or locked */
int free_twosided_helpmutex_orig(int inum_orig)
{
    int mproc=-NProcs_Work();
    int number=1;

    if (DEBUG_) printf("%5d: In free_twosided_helpmutex_orig begin: original mutex num=%d\n",ProcID(),inum_orig);

    if (inum_orig <0 || inum_orig >= MAX_TWOSIDED_HELPMUTEX_ARRAYS ) {
       fprintf(stderr,"ERROR in free_twosided_helpmutex_orig: over range! original mutex num=%d, max num=%d\n",inum_orig,MAX_TWOSIDED_HELPMUTEX_ARRAYS);
       exit (1);
    }
    else if ( twosided_helpmutex_index[inum_orig].actv==0 ) {
       fprintf(stderr,"ERROR in free_twosided_helpmutex_orig:original mutex %d has not activated for use.\n",inum_orig);
       exit (1);
    }
    else if ( twosided_helpmutex_index[inum_orig].lock==1 ) {
       fprintf(stderr,"WARNING in free_twosided_helpmutex_orig: original mutex %d is still locked before freed. Now unlocking and freeing it...\n",inum_orig);
       unlock_twosided_helpmutex_orig(inum_orig);
       twosided_helpmutex_collect(mproc, inum_orig, number);
    }
    else {
       twosided_helpmutex_collect(mproc, inum_orig, number);
    }

    if (DEBUG_) printf("%5d: In free_twosided_helpmutex_orig end: original mutex=%d\n",ProcID(),inum_orig);

    return 0;
}

/* allocates a set containing the number of mutexes. Only one set of mutexes can exist at a time. Mutexes can be
allocated and freed as many times as needed. Mutexes are numbered: 0, ..., number-1. */
/* It is only called by NON-GA routines, CAN NOT be called inside mpiga operations */
int alloc_twosided_helpmutexes(int number)
{
    int mproc=NProcs_Work();
    int inum=MAX_TWOSIDED_HELPMUTEX_GA;

    if (DEBUG_) printf("%5d: In alloc_twosided_helpmutexes: begin.\n",ProcID());
    if ( number <=0 || number > MAX_TWOSIDED_HELPMUTEX_NONGA ) {
       fprintf(stderr,"ERROR in alloc_twosided_helpmutexes: over range! number=%d, max num=%d\n",number,MAX_TWOSIDED_HELPMUTEX_NONGA);
       exit (1);
    }
    if ( twosided_helpmutex_num > 0 ) {
       fprintf(stderr,"ERROR in alloc_twosided_helpmutexes: a set of mutexes have already allocated for NONGA.\n");
       exit (1);
    }

    twosided_helpmutex_collect(mproc, inum, number);

    if (DEBUG_) printf("%5d: In alloc_twosided_helpmutexes end:  %d mutexes have been allocated.\n",ProcID(),number);

    return 0;
}

/* frees a set of mutexes allocated with alloc_twosided_helpmutexes.*/
/* It is only called by NON-GA routines, CAN NOT be called inside mpiga operations */
int free_twosided_helpmutexes()
{
    int mproc=-NProcs_Work();
    int inum=MAX_TWOSIDED_HELPMUTEX_GA;
    int number;

    if (DEBUG_) printf("%5d: In free_twosided_helpmutexes begin: twosided_helpmutex_num=%d\n",ProcID(),twosided_helpmutex_num);

    if (twosided_helpmutex_num<=0) {
       fprintf(stderr,"ERROR in free_twosided_helpmutexes: no mutex is needed to be freed. N mutex=%d\n",twosided_helpmutex_num);
       exit (1);
    }
    number=twosided_helpmutex_num;
    twosided_helpmutex_collect(mproc, inum, number);

    if (DEBUG_) printf("%5d: In free_twosided_helpmutexes end: now twosided_helpmutex_num=%d\n",ProcID(),twosided_helpmutex_num);

    return 0;
}


/* Initialization for twosided_helpmutex, allocate memory for global mutex structure twosided_helpmutex_index */
/* It is called just before install_twosided_nxtval */
void initialize_twosided_helpmutexes()
{
    int i;

    if (DEBUG_) printf("%5d: In mpi_nxtval: initialize_twosided_helpmutexes begin.\n",ProcID());
    if ( twosided_helpmutex_index != NULL ) {
       fprintf(stderr,"In mpi_nxtval: initialize_twosided_helpmutexes ERROR:  global mutex structure twosided_helpmutex_index have already existed.\n");
       exit (1);
    }
    else {
       twosided_helpmutex_index=(twosided_helpmutex_array_t *)malloc(sizeof(twosided_helpmutex_array_t)*MAX_TWOSIDED_HELPMUTEX_ARRAYS);
       if(twosided_helpmutex_index==NULL){
          fprintf(stderr,"In mpi_nxtval: initialize_twosided_helpmutexes ERROR: failed to allocate memory for global mutex structure twosided_helpmutex_index.\n");
          exit (1);
       }
    }
    for(i=0;i<MAX_TWOSIDED_HELPMUTEX_ARRAYS; i++) {
       twosided_helpmutex_index[i].actv = 0;
       twosided_helpmutex_index[i].lock = 0;
       twosided_helpmutex_index[i].mutex = NULL;
    }
    twosided_helpmutex_num=0;
    if (DEBUG_) printf("%5d: In mpi_nxtval: initialize_twosided_helpmutexes end. Number of initialized mutexes = %d\n",ProcID(),MAX_TWOSIDED_HELPMUTEX_ARRAYS);
}

/* Finalization for twosided_helpmutex, free memory for global mutex structure twosided_helpmutex_index */
/* It is called just after finalize_twosided_nxtval */
void finalize_twosided_helpmutexes()
{
    int i;

    if (DEBUG_) printf("%5d: In mpi_nxtval: finalize_twosided_helpmutexes begin.\n",ProcID());

    if(!twosided_helpmutex_index){
       if (DEBUG_) printf("%5d: In mpi_nxtval: finalize_twosided_helpmutexes: global mutex structure lists twosided_helpmutex_index do not exist.\n",ProcID());
       return;
    }
    for(i=0;i<MAX_TWOSIDED_HELPMUTEX_ARRAYS; i++) {
       if ( twosided_helpmutex_index[i].lock==1 ) {
          fprintf(stderr,"In mpi_nxtval: finalize_twosided_helpmutexes ERROR: mutex %d is still locked before freed.\n",i);
          exit (1);
       }
       if ( twosided_helpmutex_index[i].mutex != NULL ) { free (twosided_helpmutex_index[i].mutex); twosided_helpmutex_index[i].mutex=NULL; }
    }
    if ( twosided_helpmutex_index != NULL ) { free (twosided_helpmutex_index); twosided_helpmutex_index=NULL; }
    twosided_helpmutex_num=0;

    if (DEBUG_) printf("%5d: In mpi_nxtval: finalize_twosided_helpmutexes end. Number of finalized mutexes = %d\n",ProcID(),MAX_TWOSIDED_HELPMUTEX_ARRAYS);
}


/* initialise twosided_helpga data structure, zeroize for pointers in TWOSIDED_HELPGA array */
void initialize_twosided_helpga()
{
    int i,size;
    twosided_helpga_data_struc=(twosided_helpga_array_t *)malloc(sizeof(twosided_helpga_array_t)*MAX_TWOSIDED_HELPGA_ARRAYS);
    if(!twosided_helpga_data_struc){
      fprintf(stderr,"In mpi_nxtval: initialize_twosided_helpga ERROR: failed to allocate memory for twosided_helpga_data_struc.\n");
      exit (1);
    }
    twosided_helpga_index = twosided_helpga_data_struc;
    for(i=0;i<MAX_TWOSIDED_HELPGA_ARRAYS; i++) {
       twosided_helpga_index[i].actv = 0;
       twosided_helpga_index[i].ptr  = NULL;
    }
    twosided_helpga_num=0;
    twosided_helpga_curmem=(long)0;
    NUMBER_OF_SERVER=TotalNumber_of_Servers();               /* get global variable NUMBER_OF_SERVER       */
    size=NUMBER_OF_SERVER;
    twosided_helpga_map=(int *)malloc(2*size*sizeof(int));        /* initialize list of lower and upper indices */
    twosided_helpga_proclist=(int *)malloc(size*sizeof(int));     /* initialize list of processes               */
    if (DEBUG_) printf("%5d: In mpi_nxtval: initialize_twosided_helpga end. Maximum helpga=%d\n",ProcID(),MAX_TWOSIDED_HELPGA_ARRAYS);
}

/* release a helpga and free the underlying memory */
int twosided_helpga_release_orig( int handle_orig )
{
    int sizeofdtype;        /* Size of MPI Datatype */
    long sizetot;
    MPIHELPGA helpga=NULL;

    if (DEBUG_) printf("%5d: In mpi_nxtval: twosided_helpga_release_orig begin. handle_orig=%d\n",ProcID(),handle_orig);

    /* check whether helpga handle_orig is out of range, and check whether it is allocated */
    if(handle_orig < 0 || handle_orig >= MAX_TWOSIDED_HELPGA_ARRAYS){
       fprintf(stderr,"%5d: In mpi_nxtval: twosided_helpga_release_orig ERROR. invalid handle_orig [%d]. Should be [ 0 -- %d ].\n",ProcID(),handle_orig,TWOSIDED_HELPGA_OFFSET);
       exit (1);
    }
    if(twosided_helpga_index[handle_orig].actv==0){
       fprintf(stderr,"%5d: In mpi_nxtval: twosided_helpga_release_orig ERROR. no memory allocated for this helpga [handle_orig=%d].\n",ProcID(),handle_orig);
       exit (1);
    }

    helpga=twosided_helpga_index[handle_orig].ptr;
    MPI_Type_size( helpga->dtype, &sizeofdtype );
    sizetot=(long)(helpga->nele)* (long)(sizeofdtype);
    if (helpga!=NULL) { /* free the memory in helpga structure */
      if ( helpga->name != NULL ) { free (helpga->name); helpga->name=NULL;}
/*      if ( helpga->ptr_buf != NULL ) { free (helpga->ptr_buf); helpga->ptr_buf=NULL;} */
      if ( helpga->ptr_buf != NULL ) { MPI_Free_mem(helpga->ptr_buf); helpga->ptr_buf=NULL;}
      if ( helpga->len != NULL ) { free (helpga->len);  helpga->len=NULL;}
      if ( helpga->len_help != NULL ) { free ( helpga->len_help);  helpga->len_help=NULL;}
      free(helpga);
      helpga=NULL;
    }
    twosided_helpga_index[handle_orig].actv=0;
    twosided_helpga_index[handle_orig].ptr=NULL;
    twosided_helpga_num--;
    twosided_helpga_curmem -= sizetot;

    if (DEBUG_) printf("%5d: In mpi_nxtval: twosided_helpga_release_orig end. handle_orig=%d\n",ProcID(),handle_orig);

    return 0;
}


/* free the memory for global twosided_helpga data structure */
void finalize_twosided_helpga()
{
    int i;
    if (DEBUG_) printf("%5d: In mpi_nxtval: finalize_twosided_helpga begin. twosided_helpga_num=%d\n",ProcID(),twosided_helpga_num);
    if (twosided_helpga_num>0) {
       for(i=0;i<MAX_TWOSIDED_HELPGA_ARRAYS; i++) {
         if ( twosided_helpga_index[i].actv == 1 )  twosided_helpga_release_orig(i);
       }
    }
    if(twosided_helpga_data_struc) free(twosided_helpga_data_struc);
    if(twosided_helpga_map) free(twosided_helpga_map);           /* free the memory for list of lower and upper indices */
    if(twosided_helpga_proclist) free(twosided_helpga_proclist); /* free the memory for list of list of processes       */
    if (DEBUG_) printf("%5d: In mpi_nxtval: finalize_twosided_helpga end.\n",ProcID());
}

/* initialization for nxtval -- called in mpiga_initialize */
void install_twosided_nxtval()
{
    int numprocs, myid;

    if (use_helper_server) {
       initialize_twosided_helpga(); /* initialise twosided_helpga data structure  */
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (use_helper_server) {
       /* in this mode all servers are hidden from the application */
       if(SR_parallel && myid == Server_of_Rank(myid)) {
         if (DEBUG_) {
           printf("%5d: In mpi_nxtval [install_twosided_nxtval]: excluding one/more processes for data servers.\n",ProcID());
           fflush(stdout);
         }
         DataHelperServer();
         if (DEBUG_) printf("%5d: In mpi_nxtval [install_twosided_nxtval]: DataHelperServer has been terminated. Call mpiga_cleanup_finalize()\n",ProcID());
         mpiga_cleanup_finalize();
         exit(EXIT_SUCCESS);
       }
    }
    else {
       if (DEBUG_) {
         printf("%5d: In mpi_nxtval [install_twosided_nxtval]: data server is disabled.\n",ProcID());
         printf("%5d: All functions relying on data server will be unavailable!\n",ProcID());
       }
    }
}

void finalize_twosided_nxtval()
{
    if (use_helper_server) {
       finalize_twosided_helpga(); /* free global twosided_helpga data structure */
    }
}


#else

void mpi_nxtval_dummy () {}

#endif
