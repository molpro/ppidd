/* ====================================================================== *\
 *          Parallel Programming Interface for Distributed Data           *
 *          ===================================================           *
 *  Copyright (C) 2008-2010 by the authors of PPIDD. All rights reserved. *
 * ---------------------------------------------------------------------- *
 * Common Fortran and C interfaces.                                       *
 * For Fortran interface, the subroutines in this file named PPIDD_XXXXX  *
 * are converted to the proper FORTRAN external by the FORT_Extern macro  *
 * and the definitions in the ppidd_fortran.h header file.                *
 * For C interface, the subroutines in this file named PPIDD_XXXXX can be *
 * only called by C program directly. Any calling by Fortran progam       *
 * should refer to the routines in the ppidd_fortran.h file.              *
 *                                                                        *
 * Written by: Manhui Wang                                                *
 * Date:       09/10/2009                                                 *
\* ====================================================================== */

/*! \file
 * \brief This file contains the code shared by both Fortran and C interface routines to PPIDD */

#ifdef PPIDD_LANG
/* PPIDD_LANG
 * 1 = C
 * 2 = Fortran
 */
#else
#error "PPIDD_LANG not defined"
#endif

#if defined(MPI2) || defined(GA_TOOLS)

#ifdef MPI2
 #include "mpiga_base.h"         /* include mpi.h, ppidd_machines.h, and ppidd_dtype.h */
 #include "mpi_utils.h"
 #include "mpi_nxtval.h"
#endif

#ifdef GA_TOOLS
 #include "ppidd_machines.h"
 #include <ga.h>
 #include <macdecls.h>

  #ifdef GA_VERSION_GE_5
    #if (GA_VERSION_MAJOR == 5) && (GA_VERSION_MINOR == 0)
      extern int ga_type_f2c(int type);
    #else
      #include <ga-papi.h>
      #define ga_type_f2c pnga_type_f2c
    #endif
  #endif

  #ifdef GA_MPI
    #include <mpi.h>
    #include "mpi_utils.h"
  #else
    #ifdef GA_TCGMSGMPI
      #include <mpi.h>
      #include "mpi_utils.h"
    #endif
    #ifdef GA_VERSION_GE_5
      #include <tcgmsg.h>
    #else
      #include <sndrcv.h>
    #endif
    #if (PPIDD_LANG == 2)
/*! \cond */
    void PBEGINF_();  /*! use fortran interface routine to avoid handling the arguments additionally */
/*! \endcond */
    #endif
  #endif
#endif
 static int MPIGA_Debug=0;

/*! For serial case */
#else
 #include "ppidd_machines.h"
#endif

#if (PPIDD_LANG == 1)
#include "ppidd_undefdtype.h"
#include "ppidd_c.h"
#endif

#if (PPIDD_LANG == 2)
#include "ppidd_dtype.h"
#include "ppidd_fortran.h"
#endif

/*! Initialize the PPIDD parallel environment
 *
 *  - For \b GA, includes initialization of TCGMSG/MPI and GA.
 *  - For \b MPI2, calls MPI_Init.
 */
#if (PPIDD_LANG == 1)
   void PPIDD_Initialize(int argc, char **argv) {
#endif
#if (PPIDD_LANG == 2)
   void ppidd_initialize_f2c(int nargs, char **args) {
#endif
#if defined(MPI2) || defined(GA_MPI)

  #ifdef MPI2
      int mpierr;
  #endif
#if (PPIDD_LANG == 2)
      int i,j,argc=nargs;
      char **argv = NULL;
      char arg[256];

   /* ========================== *\
      Get command line arguments
   \* ========================== */
      if(argc>0) {
         argv = malloc(argc*sizeof(char*));
         for(i=0; i<argc; ++i) {
           argv[i] = (char *) strdup(args[i]);
           if(MPIGA_Debug)printf("In ppidd_initialize: argc=%d, argv[%d]=%s\n",argc,i,argv[i]);
         }
      }
      else {
/* If no argument is provided, use the dummy command arguments, since MPI_Init can't work if arguments are null for MPICH1.
   For MPICH2 and other MPI2 library, no such limitation */
         argc=1;
         i=0;
         argv = malloc(argc*sizeof(char*));
         strcpy(arg,"DummyCommand");
         argv[i] = (char *) strdup(arg);
         if(MPIGA_Debug)printf("In ppidd_initialize: USE DUMMY ARGUMENTS. argc=%d, argv[%d]=%s\n",argc,i,argv[i]);
      }
#endif

   /* ================ *\
      Initialize MPIGA
   \* ================ */
  #ifdef MPI2
      mpierr=mpiga_initialize(&argc,&argv);
      mpi_test_status("PPIDD_Initialize:",mpierr);
  #endif
  #ifdef GA_MPI
      MPI_Init(&argc, &argv);                       /* initialize MPI */
  #endif
#endif
#ifdef GA_TOOLS
  #if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
/* call c interface routine PBEGIN_(argc, argv) since it will get the right arguments from c main program */
/* call fortran interface routine PBEGINF_() instead of c interface routine PBEGIN_(argc, argv) to avoid handling the arguments additionally */
    #ifdef GA_VERSION_GE_5
     #if (PPIDD_LANG == 1)
      tcg_pbegin(argc, argv);                    /* initialize TCGMSG */
     #endif
     #if (PPIDD_LANG == 2)
      tcg_pbeginf();                             /* initialize TCGMSG */
     #endif
    #else
     #if (PPIDD_LANG == 1)
      PBEGIN_(argc, argv);                       /* initialize TCGMSG */
     #endif
     #if (PPIDD_LANG == 2)
      PBEGINF_();                                /* initialize TCGMSG */
     #endif
    #endif
  #endif
      GA_Initialize();                            /* initialize GA */
#endif

#if (PPIDD_LANG == 2)
   /* ================================================================== *\
      Free the memory for argv[i] and argv itself after having used them
   \* ================================================================== */
   /* strdup implies copying a string to a particular location that will be automatically created by the function.
    * The allocated memory for argv[i] needs to be explicitly freed before leaving this function */
#if defined(MPI2) || defined(GA_MPI)
      for (i = 0; i < argc; i++) {
         if (argv[i]) free(argv[i]);
      }
   /* Also need to free argv itself */
      if (argv) free(argv);
  #ifdef MPI2
      if(MPIGA_Debug)printf("%5d: In ppidd_share: end of ppidd_initialize.\n",ProcID());
  #endif
#endif
#endif
   }

#if (PPIDD_LANG == 2)
   void PPIDD_Initialize(void) {
    int i=0;
    char **a = NULL;
    ppidd_initialize_f2c(i,a);
   }
#endif

#if (PPIDD_LANG == 1)
/* Both fortint and fortlogical are defined as int. */
/* FALSE and TRUE are (int)0 and (int)1, respectively. */
#endif


/*! Initialize the PPIDD data structure
 *
 *  - For \b GA, does nothing
 *  - For \b MPI2, Initialize global data structure and set helper server.
 */
   void PPIDD_Initialize_data(void) {
#ifdef MPI2
      mpiga_initialize_data();
      if(MPIGA_Debug)printf("%5d: In ppidd_share: end of ppidd_initialize_data.\n",ProcID());
#endif
   }


/*! Return communicator for worker process group, which excludes any helper process.
 * For simplicity, the returned data type is set as integer. So for C calling, it may have to be converted
 * into MPI_Comm by MPI_Comm_f2c afterward if needed.
 *
 *  - For \b MPI2, Return communicator for worker process group.
 *  - For \b GA and serial cases, should not be called.
 */
   fortint PPIDD_Worker_comm(void) {
#ifdef MPI2
      MPI_Comm mycomm=mpigv(Compute_comm);
      MPI_Fint fcomm;

/* test whether worker communicator contains all the processes, if so then return MPI_COMM_WORLD */
      int np_all, np_worker=mpigv(nprocs);
      MPI_Comm_size(MPI_COMM_WORLD, &np_all);
      if(np_all==np_worker) mycomm=MPI_COMM_WORLD;

      fcomm=MPI_Comm_c2f(mycomm);
      return (fortint)fcomm;
#elif defined(GA_TCGMSGMPI) || defined(GA_MPI)
      MPI_Fint fcomm;
      fcomm=MPI_Comm_c2f(MPI_COMM_WORLD);
      return (fortint)fcomm;
#else
      fprintf(stderr," ERROR: PPIDD_Worker_comm should not be called in NON-MPI2 cases.\n");
      return (fortint)-1;
#endif
   }


/*! Terminate the PPIDD parallel environment.
 *
 *  - For \b GA, tidy up global arrays and TCGMSG/MPI, analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_terminate
 *  - For \b MPI2, tidy up some associated resources and call MPI_Finalize.
 */
   void PPIDD_Finalize(void) {
#ifdef MPI2
      mpiga_terminate();
#endif
#ifdef GA_TOOLS
      GA_Terminate();
#ifdef GA_MPI
      MPI_Finalize();
#endif
#if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
    #ifndef GA_VERSION_GE_5
      PEND_();
    #else
      tcg_pend();
    #endif
#endif
#endif
   }


/*! \cond */
/*  Detect whether MA is used for allocation of GA memory.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_uses_ma
 *  - \b MPI2 always returns <tt>.false.</tt>
 */
   void PPIDD_Uses_ma(fortlogical *ok) {
#ifdef MPI2
      *ok = FALSE;
#endif
#ifdef GA_TOOLS
      if(MPIGA_Debug)printf("In ppidd_uses_ma: sizeof(fortlogical)=%d, sizeof(fortint)=%d, sizeof(Integer)=%d\n",
                             (int)sizeof(fortlogical),(int)sizeof(fortint),(int)sizeof(Integer));
      if(MPIGA_Debug)printf("In ppidd_uses_ma: sizeof(double)=%d,  sizeof(DoublePrecision)=%d\n",(int)sizeof(double),(int)sizeof(DoublePrecision));

      if (GA_Uses_ma()) *ok = TRUE;
      else *ok = FALSE;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *ok = FALSE;
#endif
  }


/*  Initialize the memory allocator.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/parsoft/ma/MAapi.html
 *  - \b MPI2 always returns <tt>.true.</tt>
 */
   void PPIDD_MA_init(fortint *dtype, fortint *stack, fortint *heap, fortlogical *ok) {
#ifdef MPI2
      *ok = TRUE;
#endif
#ifdef GA_TOOLS
      int istack=(int)*stack;
      int iheap=(int)*heap;
      int gadtype=-1;
      char *errmsg;

      switch((int)*dtype){
      case 0:
              gadtype=MT_F_INT;
              break;
      case 1:
              gadtype=MT_F_DBL;
              break;
      default:
              errmsg=strdup(" In PPIDD_MA_Init: wrong data type ");
              GA_Error(errmsg,(int)*dtype);
              free(errmsg);
      }
      if( MA_init(gadtype, istack, iheap)) *ok = TRUE;
      else *ok = FALSE;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *ok = TRUE;
#endif
  }
/*! \endcond */

/*! Return an elapsed time on the calling process.
 *
 * Returns a floating-point number of seconds, representing elapsed wall-clock time since an arbitrary time in the past.
 * The times returned are local to the node that called them. There is no requirement that different nodes return the same time.
 * This is a local operation.
 *
 *  - \b GA calls GA_Wtime, which is available only in release 4.1 or greater, http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_wtime
 *  - \b MPI2 calls MPI_Wtime
 */
   void PPIDD_Wtime(double *ctime) {
#ifdef MPI2
      *ctime = MPI_Wtime();
#endif
#ifdef GA_TOOLS
      *ctime = GA_Wtime();
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *ctime = (double)0;
#endif
  }


/*! Print an error message and abort the program execution.
 *
 *  - \b GA calls GA_Error, http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_error
 *  - For \b MPI2, prints error, and then calls MPI_Abort.
 */
   void PPIDD_Error(char *message
/*! \cond */
#if defined(FORTCL_NEXT)
	       ,fortintc lx
#endif
/*! \endcond */
	       ,fortint *code
/*! \cond */
#if defined(FORTCL_END)
	       ,fortintc lx
#endif
/*! \endcond */
      ) {
      char *msg2, *p;
      int icode=(int)*code;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(message);
#endif
      strncpy((msg2=(char *)calloc(lxi+1,1)),message,lxi);
      for (p=msg2+lxi;p>=msg2;p--) if (*p==' ')*p=(char)0;

#ifdef MPI2
      MPIGA_Error(msg2, icode);
#endif
#ifdef GA_TOOLS
      GA_Error(msg2, icode);
#endif
#if defined(MPI2) || defined(GA_TOOLS)
      free(msg2);
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      fprintf(stdout," %s %d (%#x).\n", msg2,icode,icode);
      fflush(stdout);
      fprintf(stderr," %s %d (%#x).\n", msg2,icode,icode);

      printf(" PPIDD_Error: now exiting...\n");
      free(msg2);
      exit(1);
#endif
   }


/*! Set the flag for data helper server and set the number of how many processes own one helper server
 *
 *  Set the helper_server flag: 1 (use); 0 (don't use). */
   void PPIDD_Helper_server(fortlogical *flag, fortint *numprocs_per_server) {
#ifdef MPI2
      if ((int)*flag) {                              /* mutilple helper servers, node helper server, and single helper server */
         use_helper_server=1;
         if ( (int)*numprocs_per_server > 1 ) {      /* reasonable mutilple helper servers, and single helper server */
            NPROCS_PER_HELPER=(int)*numprocs_per_server;
         }
         else if ( (int)*numprocs_per_server == 0 ) {/* node helper server: one helper server on every node */
           int mpinp;
           MPI_Comm_size(MPI_COMM_WORLD, &mpinp);
           if (NNODES_SYMMETRY) {
             NPROCS_PER_HELPER=mpinp/NUM_TOTAL_NNODES;
             if (NPROCS_PER_HELPER==1) {             /* node helper server: if NPROCS_PER_HELPER==1, then use only one single helper server */
               if(NUM_TOTAL_NNODES>1)fprintf(stdout,"%5d: WARNING: only one process on each node. Will use only one single helper server for all processes!\n", ProcID());
               NPROCS_PER_HELPER=99999999;
             }
           }
           else {                                    /* node helper server: if not all the nodes are symmetric, then use only one single helper server */
             fprintf(stderr,"%5d: WARNING: not all the nodes are symmetric. Will use only one single helper server for all processes!\n", ProcID());
             NPROCS_PER_HELPER=99999999;
           }
         }
         else {                                      /* unreasonable mutilple helper servers, then only use single helper server */
            fprintf(stderr,"%5d: WARNING: nprocs_per_server=%d is unreasonable. Will use only single helper server for all processes!\n", ProcID(),(int)*numprocs_per_server);
            NPROCS_PER_HELPER=99999999;
         }
      }
      else {                                         /* no helper server */
         use_helper_server=0;
         NPROCS_PER_HELPER=-1;
      }
#endif
   }


/*! Get the next shared counter number(helper process should be enabled).
 *
 *  Increment a counter by 1 and returns the counter value (0, 1, ...). */
   void PPIDD_Nxtval(fortint *numproc, fortint *val) {
#ifdef MPI2
      if (use_helper_server==0) {
        fprintf(stderr,"%5d: ERROR: Attemp to call NXTVAL routine without helper process!\n", ProcID());
        MPI_Abort(mpigv(Compute_comm),911);
      }
      else {
        int mproc = (int) *numproc;
        *val= (fortint) NXTVAL(&mproc);
      }
#endif
#ifdef GA_TOOLS
  #if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
      long mproc = (long) *numproc;
    #ifndef GA_VERSION_GE_5
      *val= (fortint) NXTVAL_(&mproc);
    #else
      *val= (fortint) tcg_nxtval(mproc);
    #endif
  #endif
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Nxtval should not be called in serial case.\n");
      exit(1);
#endif
   }



/*!  Determine the total number of processes available (including helper process if there is one).
 *
 *  - \b GA calls GA_Nnodes,
 *  http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_nnodes
 *  - \b MPI2 calls MPI_Comm_size for communicator MPI_COMM_WORLD.
 */
   void PPIDD_Size_all(fortint *np) {
#ifdef MPI2
      int mpinp;
      MPI_Comm	mpicomm=MPI_COMM_WORLD;

      MPI_Comm_size(mpicomm, &mpinp);
      *np = (fortint) mpinp;
#endif
#ifdef GA_TOOLS
      *np = (fortint)GA_Nnodes();
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *np = (fortint)1;
#endif
   }


/*  Determine the number of compute processes.
 *
 *  - \b GA calls GA_Nnodes,
 *  http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_nnodes
 *  - \b MPI2 calls MPI_Comm_size for computational communicator.
 */
   void PPIDD_Size(fortint *np) {
#ifdef MPI2
      int mpinp;
      MPI_Comm	mpicomm=mpigv(Compute_comm);

      MPI_Comm_size(mpicomm, &mpinp);
      *np = (fortint) mpinp;
#endif
#ifdef GA_TOOLS
      *np = (fortint)GA_Nnodes();
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *np = (fortint)1;
#endif
   }


/*! Determine the rank of the calling process.
 *
 *  - \b GA calls GA_Nodeid,
 *  http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_nodeid
 *  - \b MPI2 calls MPI_Comm_rank in computational communicator.
 */
   void PPIDD_Rank(fortint *me) {
#ifdef MPI2
      int mpime;
      MPI_Comm	mpicomm=mpigv(Compute_comm);

      MPI_Comm_rank(mpicomm, &mpime);
      *me = (fortint) mpime;
#endif
#ifdef GA_TOOLS
      *me = (fortint)GA_Nodeid();
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *me = (fortint)0;
#endif
   }


/*! Initialize tracing of completion status of data movement operations.
 *
 *  - \b GA calls GA_Init_fence, http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_init_fence
 *  - \b MPI2 does nothing
 */
   void PPIDD_Init_fence(void) {
#ifdef GA_TOOLS
      GA_Init_fence();
#endif
   }


/*! Block the calling process until all data transfers complete.
 *
 *  - \b GA calls GA_Fence, http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_fence
 *  - \b MPI2 does nothing
 */
   void PPIDD_Fence(void) {
#ifdef GA_TOOLS
      GA_Fence();
#endif
   }



#if defined(MPI2) || defined(GA_MPI)
/* =================== nonblocking mpi message list ================= */
#define MAX_MPIQ_LEN 128      /* Maximum no. of outstanding messages */
/* In application programs, nonblocking send/recv may be never used, so here set it to a small number. It could be increased if necessary. */
/*! \cond */
static struct msg_mpiq_struct{
  MPI_Request request;
  long         node;
  long         type;
  long         lenbuf;
  long         snd;
  long         from;
} msg_mpiq[MAX_MPIQ_LEN];
/*! \endcond */

static int n_in_msg_mpiq=0;
/* =================================================================== */
#endif


/*! Blocking/nonblocking send.
 *
 *  - \b GA calls SND_.
 *  - \b MPI2 calls MPI_Send ( sync is 1) or MPI_Isend ( sync is 0).
 */
   void PPIDD_Send(void *buf,fortint *count,fortint *dtype,fortint *dest,fortint *sync) {
#if defined(MPI2) || defined(GA_MPI)
  #ifdef MPI2
      MPI_Comm	mpicomm=mpigv(Compute_comm);
  #endif
  #ifdef GA_MPI
      MPI_Comm	mpicomm = MPI_COMM_WORLD;
  #endif
      int mpicount=(int)*count;
      int mpidest=(int)*dest;
      int mpitag=(int)*dtype;
      int mpisync=(int)*sync;
      int mpilenbuf,mpierr;
      MPI_Datatype mpidtype;
      int sizempidtype;

      mpiga_type_f2cmpi((int)*dtype,&mpidtype,&sizempidtype);
      mpilenbuf=mpicount*sizempidtype;

      if (MPIGA_Debug) {
         printf("PPIDD_SEND: node %d sending to %d, len(bytes)=%d, mes tag=%d, sync=%d\n",
                 ProcID(), mpidest, mpilenbuf, mpitag, mpisync);
         fflush(stdout);
      }

      if (mpisync) {
         mpierr=MPI_Send(buf,mpilenbuf,MPI_CHAR,mpidest,mpitag,mpicomm);
         mpi_test_status("PPIDD_SEND: SEND:",mpierr);
      }
      else {
         if (n_in_msg_mpiq >= MAX_MPIQ_LEN) {
            MPIGA_Error("PPIDD_SEND: nonblocking SEND: overflowing async Queue limit",n_in_msg_mpiq);
	 }
         mpierr = MPI_Isend(buf, mpilenbuf, MPI_CHAR,mpidest, mpitag,mpicomm,
                     &msg_mpiq[n_in_msg_mpiq].request);
         mpi_test_status("PPIDD_SEND: nonblocking SEND:",mpierr);

         msg_mpiq[n_in_msg_mpiq].node   =(long) mpidest;
         msg_mpiq[n_in_msg_mpiq].type   =(long) *dtype;
         msg_mpiq[n_in_msg_mpiq].lenbuf =(long) mpilenbuf;
         msg_mpiq[n_in_msg_mpiq].snd = (long)1;
      }
#endif
#if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
      long gadtype=(long)-1;
      long galenbuf;
      long gadest=(long)*dest;
      long gasync=(long)*sync;
      size_t ctype=8;
      char *errmsg;

      switch((int)*dtype){
      case 0:
              gadtype=MT_F_INT;
              ctype=sizeof(fortint);
              break;
      case 1:
              gadtype=MT_F_DBL;
              ctype=sizeof(double);
              break;
      default:
              errmsg=strdup(" In PPIDD_Send: wrong data type ");
              GA_Error(errmsg,(int)*dtype);
              free(errmsg);
      }
      galenbuf=(*count)*ctype;
    #ifndef GA_VERSION_GE_5
      SND_(&gadtype, buf, &galenbuf, &gadest, &gasync);
    #else
      tcg_snd(gadtype, buf, galenbuf, gadest, gasync);
    #endif
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Send should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Blocking/nonblocking receive.
 *
 *  - \b GA calls RCV_.
 *  - \b MPI2 calls MPI_Recv ( sync is 1) or MPI_Irecv( sync is 0).
 */
   void PPIDD_Recv(void *buf,fortint *count,fortint *dtype,fortint *source,fortint *lenreal,fortint *sourcereal,fortint *sync) {
#if defined(MPI2) || defined(GA_MPI)
  #ifdef MPI2
      MPI_Comm	mpicomm=mpigv(Compute_comm);
  #endif
  #ifdef GA_MPI
      MPI_Comm	mpicomm = MPI_COMM_WORLD;
  #endif
      int mpicount=(int)*count;
      int mpitag=(int)*dtype;
      int mpisource=(int)*source;
      int mpisync=(int)*sync;
      int mpinode,mpilenbuf,mpierr;
      MPI_Status status;
      MPI_Request request;
      MPI_Datatype mpidtype;
      int sizempidtype;

      mpiga_type_f2cmpi((int)*dtype,&mpidtype,&sizempidtype);
      mpilenbuf=mpicount*sizempidtype;

      if (mpisource == -1)
         mpinode = MPI_ANY_SOURCE;
      else
         mpinode = mpisource;

      if (MPIGA_Debug) {
         printf("PPIDD_Recv: node %d receving from %d, len(bytes)=%d, mes tag=%d, sync=%d\n",
                 ProcID(), mpisource, mpilenbuf, mpitag, mpisync);
         fflush(stdout);
      }

      if(mpisync==0){
         if (n_in_msg_mpiq >= MAX_MPIQ_LEN) {
            MPIGA_Error("PPIDD_Recv: nonblocking RECV: overflowing async Queue limit", n_in_msg_mpiq);
	 }
         mpierr = MPI_Irecv(buf,mpilenbuf,MPI_CHAR,mpinode,mpitag,mpicomm,&request);
         mpi_test_status("PPIDD_Recv: nonblocking RECV:",mpierr);

         *sourcereal = (fortint) mpinode;          /* Get source node  */
         *lenreal = (fortint) (-1);
         msg_mpiq[n_in_msg_mpiq].request = request;
         msg_mpiq[n_in_msg_mpiq].node   = (long)*source;
         msg_mpiq[n_in_msg_mpiq].type   = (long)*dtype;
         msg_mpiq[n_in_msg_mpiq].lenbuf = (long)mpilenbuf;
         msg_mpiq[n_in_msg_mpiq].snd = (long)0;
         n_in_msg_mpiq++;
      }
      else{
         mpierr = MPI_Recv(buf,mpilenbuf,MPI_CHAR,mpinode,mpitag,mpicomm,&status);
         mpi_test_status("PPIDD_RECV: RECV:",mpierr);
         mpierr = MPI_Get_count(&status, MPI_CHAR, &mpilenbuf);
         mpi_test_status("PPIDD_RECV: Get_count:",mpierr);
         *sourcereal = (fortint)status.MPI_SOURCE;
         *lenreal    = (fortint)mpilenbuf;
      }
#endif
#if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
      long gadtype=(long)-1;
      long galenbuf;
      long gasource=(long)*source;
      long gasync=(long)*sync;
      long galenreal;
      long gasourcereal;
      size_t ctype=8;
      char *errmsg;

      switch((int)*dtype){
      case 0:
              gadtype=MT_F_INT;
              ctype=sizeof(fortint);
              break;
      case 1:
              gadtype=MT_F_DBL;
              ctype=sizeof(double);
              break;
      default:
              errmsg=strdup(" In PPIDD_Recv: wrong data type ");
              GA_Error(errmsg,(int)*dtype);
              free(errmsg);
      }

      galenbuf=(*count)*ctype;
    #ifndef GA_VERSION_GE_5
      RCV_(&gadtype, buf, &galenbuf, &galenreal, &gasource, &gasourcereal, &gasync);
    #else
      tcg_rcv(gadtype, buf, galenbuf, &galenreal, gasource, &gasourcereal, gasync);
    #endif
      *lenreal    = (fortint)galenreal;
      *sourcereal = (fortint)gasourcereal;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Recv should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Wait for completion of all asynchronous send/receive.
 *
 * ignores nodesel !! */
/*!
 *  - \b GA (with TCGMSG) analogous to WAITCOM
 *  - \b MPI2 calls MPI_Wait for all asynchronous requests.
 */
   void PPIDD_Wait(fortint *nodesel) {
#if defined(MPI2) || defined(GA_MPI)
      int mpierr,i;
      MPI_Status status;

      for (i=0; i<n_in_msg_mpiq; i++){
         if (MPIGA_Debug) {
            printf("PPIDD_Wait: node %d waiting for msg to/from %ld, #%d\n", ProcID(), msg_mpiq[i].node, i);
            fflush(stdout);
         }
         mpierr = MPI_Wait(&msg_mpiq[i].request, &status);
         mpi_test_status("PPIDD_Wait:",mpierr);
      }
      n_in_msg_mpiq = 0;
#endif
#if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
      long  ganode =  (long) *nodesel;

    #ifndef GA_VERSION_GE_5
      WAITCOM_(&ganode);
    #else
      tcg_waitcom(ganode);
    #endif
#endif
   }


/*! Detect whether a message of the specified type is available.
 *
 *  Return <tt>.true.</tt> if the message is available, otherwise <tt>.false.</tt>.
 *
 *  - \b GA (with TCGMSG) analogous to PROBE
 *  - \b MPI2 calls MPI_Iprobe
 */
   void PPIDD_Iprobe(fortint *tag,fortint *source,fortlogical *ok) {
#if defined(MPI2) || defined(GA_MPI)
  #ifdef MPI2
      MPI_Comm	mpicomm=mpigv(Compute_comm);
  #endif
  #ifdef GA_MPI
      MPI_Comm	mpicomm = MPI_COMM_WORLD;
  #endif
      int mpitag=(int)*tag;
      int mpisource,mpierr,flag;
      MPI_Status status;

      mpisource = (*source < 0) ? MPI_ANY_SOURCE  : (int) *source;
      mpierr = MPI_Iprobe(mpisource, mpitag, mpicomm, &flag, &status);
      mpi_test_status("PPIDD_Iprobe:",mpierr);
      if(flag) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#if defined(GA_TCGMSG) || defined(GA_TCGMSG5) || defined(GA_TCGMSGMPI)
      long gatag=(long)*tag;
      long gasource=(long)*source;
      long gaflag;
    #ifndef GA_VERSION_GE_5
      gaflag=PROBE_(&gatag, &gasource);
    #else
      gaflag=tcg_probe(gatag, gasource);
    #endif
      if(gaflag) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
    *ok = FALSE ;
#endif
   }


/*! Broadcast a message from the root process to all other processes.
 *
 *  Collective operation.
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_brdcst
 *  - \b MPI2 calls MPI_Bcast
 *
 *  - \c type=0 : Fortran integer and logical types
 *  - \c type=1 : Fortran double precision type */
   void PPIDD_BCast(void *buffer,fortint *count,fortint *type,fortint *root) {
#ifdef MPI2
      int mpicount=(int)*count;
      int dtype =(int)*type;
      int mpiroot=(int)*root;
      MPI_Datatype mpidtype;
      int sizempidtype;
      MPI_Comm mpicomm;
      int mpierr;

      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      mpicomm=mpigv(Compute_comm);

      mpierr=MPI_Bcast(buffer,mpicount,mpidtype,mpiroot,mpicomm);
      mpi_test_status("PPIDD_BCast:",mpierr);
#endif
#ifdef GA_TOOLS
      int gacount=(int)*count;
      int dtype =(int)*type;
      int garoot=(int)*root;
      int galenbuf;  /* in bytes */
      size_t ctype=8;
      char *errmsg;

      switch(dtype){
      case 0:
              ctype=sizeof(fortint);
              break;
      case 1:
              ctype=sizeof(double);
              break;
     default:
              errmsg=strdup(" In PPIDD_BCast: wrong data type ");
              GA_Error(errmsg,dtype);
              free(errmsg);
      }
      galenbuf=ctype*gacount;
      GA_Brdcst(buffer, galenbuf, garoot);
#endif
   }


/*! Synchronize processes and ensure all have reached this routine.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_sync
 *  - \b MPI2 calls MPI_Barrier
 */
   void PPIDD_Barrier(void) {
#ifdef MPI2
      int mpierr;

      mpierr=MPI_Barrier(mpigv(Compute_comm));
      mpi_test_status("PPIDD_Barrier:",mpierr);
#endif
#ifdef GA_TOOLS
      GA_Sync();
#endif
   }



/*! Combine values from all processes and distribute the result back to all processes.
 *
 *  - \b GA analogous to GA_dgop and GA_igop, http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_dgop
 *  - For \b MPI2, it is realized by MPI_Allreduce
 *
 *  - \c type=0 : Fortran Integer
 *  - \c type=1 : Fortran Double Precision */
   void PPIDD_Gsum(fortint *type,void *buffer,fortint *len, char *op
/*! \cond */
#if defined(FORTCL_NEXT)
	       ,fortintc lx
#endif
#if defined(FORTCL_END)
	       ,fortintc lx
#endif
/*! \endcond */
      ) {
#if defined(MPI2) || defined(GA_TOOLS)
      int dtype=(int)*type;
#endif
      char *op2, *p;
      int lxi;
#ifdef MPI2
      MPI_Datatype mpidtype;
      int sizempidtype;
      int mpilen=(int)*len;
#endif
#ifdef GA_TOOLS
      int buflen=(int)*len;
      int gadtype_f=-1,gadtype_c;
      char *errmsg;
#endif

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(op);
#endif
      strncpy(op2=(char *)calloc(lxi+1,1),op,lxi);
      for (p=op2+lxi;p>=op2;p--) if (*p==' ')*p=(char)0;

#ifdef MPI2
      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      MPI_GSum(mpidtype,buffer,mpilen, op2);
#endif
#ifdef GA_TOOLS
      switch(dtype){
      case 0:
              gadtype_f=MT_F_INT;
              break;
      case 1:
              gadtype_f=MT_F_DBL;
              break;
      default:
              errmsg=strdup(" In PPIDD_Gsum: wrong data type ");
              GA_Error(errmsg,dtype);
              free(errmsg);
      }
      gadtype_c=ga_type_f2c(gadtype_f);
      GA_Gop(gadtype_c, buffer, buflen, op2);
#endif
      free(op2);
   }




/*! Create an array by following the user-specified distribution and return integer handle representing the array.
 *
 * Irregular distributed array data are stored across the distributed processes.
 *  - \c datatype=0 : Fortran integer and logical types
 *  - \c datatype=1 : Fortran double precision type

 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_create_irreg
 */
   void PPIDD_Create_irreg(char *name
/*! \cond */
#if defined(FORTCL_NEXT)
	       ,fortintc lx
#endif
/*! \endcond */
	       ,fortint *lenin, fortint *nchunk, fortint *datatype, fortint *storetype, fortint *handle, fortlogical *ok
/*! \cond */
#if defined(FORTCL_END)
	       ,fortintc lx
#endif
/*! \endcond */
      ) {
#ifdef MPI2
      int *mpilenin,i,mpierr;
      int mpinchunk=(int)*nchunk;
      int dtype=(int)*datatype;
      int stype=(int)*storetype;
      MPI_Datatype mpidtype;
      int sizempidtype;
      int mpihandle;
      char *name2;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
      strncpy((name2=(char *)malloc(lxi+1)),name,lxi);
      name2[lxi]=(char)0;
      for(i=lxi-1; (i>=0 && name2[i]==' '); i--) name2[i]=(char)0;
      mpilenin=(int *)malloc(mpinchunk*sizeof(int));
      for (i=0;i<mpinchunk;i++) mpilenin[i]=(int)lenin[i];
      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      if (use_helper_server==0) {
        mpierr=mpiga_create_irreg(name2, mpilenin, mpinchunk, mpidtype, &mpihandle);
      }
      else {
        if (stype==0)
          mpierr=mpiga_create_irreg(name2, mpilenin, mpinchunk, mpidtype, &mpihandle);
        else {
          int mproc=0;
          mpierr=twosided_helpga_create_irreg(mproc, mpilenin, mpinchunk, &mpihandle, name2, mpidtype);
        }
      }
      free(name2);
      free(mpilenin);
      *handle=(fortint)mpihandle;
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int dtype=(int)*datatype;
      int gadtype=-1;
      int ndim=1;
      int nblock=(int)*nchunk;
      int *dims, *block, *map;
      int np;
      int i,iad,totlen;
      int gahandle;
      char *name2;
      char *errmsg;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
      strncpy((name2=(char *)malloc(lxi+1)),name,lxi);
      name2[lxi]=(char)0;
      for(i=lxi-1; (i>=0 && name2[i]==' '); i--) name2[i]=(char)0;

      switch(dtype){
      case 0:
              gadtype=MT_F_INT;
              break;
      case 1:
              gadtype=MT_F_DBL;
              break;
      default:
              errmsg=strdup(" In PPIDD_Create_irreg: wrong data type ");
              GA_Error(errmsg,dtype);
              free(errmsg);
      }

      dims=(int *)malloc(ndim*sizeof(int));
      block=(int *)malloc(ndim*sizeof(int));
      np = GA_Nnodes();
/* map[np] or map[nblock] ? */
      map=(int *)malloc(np*sizeof(int));

      for(i=0;i<ndim;i++) block[i]=nblock;

      for(iad=0,i=0;i<nblock;i++){
        map[i]=iad;
        iad=iad+(int)lenin[i];
      }
      for(i=nblock;i<np;i++) map[i]=iad;

      totlen=iad;
      for(i=0;i<ndim;i++) dims[i]=totlen;

      gahandle=NGA_Create_irreg(gadtype, ndim, dims, name2, block, map);

      free(name2);
      free(dims);
      free(block);
      free(map);

      *handle=(fortint)gahandle;
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Create_irreg should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Create an array using the regular distribution model and return integer handle representing the array.
 *
 *  - \c datatype=0  : Fortran integer and logical types
 *  - \c datatype=1  : Fortran double precision type
 *  - \c storetype=0 : Normal distributed array stored across the distributed processes
 *  - \c storetype>=1: Low-latency array stored on one or more helper processes (effective only when helper process is enabled). \c storetype is advisory: the underlying implementation will use up to \c storetype helpers.

 *  - For \b GA, storetype doesn't take effect, and data are always stored across the distributed processes, analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_create
 *  - For \b MPI2, the library can presently be built with zero or one (default) helpers.
 *       When helper process is disabled, \c storetype doesn't take effect, and data are always stored across the distributed processes.
 */
   void PPIDD_Create(char *name
/*! \cond */
#if defined(FORTCL_NEXT)
	       ,fortintc lx
#endif
/*! \endcond */
	       ,fortint *lentot, fortint *datatype, fortint *storetype, fortint *handle, fortlogical *ok
/*! \cond */
#if defined(FORTCL_END)
	       ,fortintc lx
#endif
/*! \endcond */
      ) {
#ifdef MPI2
      int mpierr;
      int mpilentot=(int)*lentot;
      int dtype=(int)*datatype;
      MPI_Datatype mpidtype;
      int sizempidtype;
      int stype=(int)*storetype;
      int mpihandle;
      char *name2, *p;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
      if(MPIGA_Debug)printf("%5d: In PPIDD_Create: sizeof(fortint)=%d,sizeof(fortintc)=%d,sizeof(fortlogical)=%d,lxi=%d\n",
	ProcID(),(int)sizeof(fortint),(int)sizeof(fortintc),(int)sizeof(fortlogical),lxi);
      strncpy((name2=(char *)calloc(lxi+1,1)),name,lxi);
      for (p=name2+lxi;p>=name2;p--) if (*p==' ')*p=(char)0;
      mpiga_type_f2cmpi(dtype,&mpidtype,&sizempidtype);
      if (use_helper_server==0) {
        mpierr=mpiga_create( name2, mpilentot, mpidtype, &mpihandle );
      }
      else {
        if (stype==0)
         mpierr=mpiga_create( name2, mpilentot, mpidtype, &mpihandle );
        else {
         int mproc=0;
         mpierr=twosided_helpga_create(mproc, mpilentot, &mpihandle, name2, mpidtype);
        }
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Create: array %s created, datatype=%d, storetype=%d\n",ProcID(),name2,stype,dtype);

      free(name2);
      *handle=(fortint)mpihandle;
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int dtype=(int)*datatype;
      int gadtype=-1;
      int ndim=1;
      int *dims, *block;
      int galentot=(int)*lentot;
      int i;
      int gahandle;
      char *name2, *p;
      char *errmsg;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
      strncpy((name2=(char *)calloc(lxi+1,1)),name,lxi);
      for (p=name2+lxi;p>=name2;p--) if (*p==' ')*p=(char)0;

      switch(dtype){
      case 0:
              gadtype=MT_F_INT;
              break;
      case 1:
              gadtype=MT_F_DBL;
              break;
      default:
              errmsg=strdup(" In PPIDD_Create: wrong data type ");
              GA_Error(errmsg,dtype);
              free(errmsg);
      }

      dims=(int *)malloc(ndim*sizeof(int));
      block=(int *)malloc(ndim*sizeof(int));

      for(i=0;i<ndim;i++) block[i]=-1;

      for(i=0;i<ndim;i++) dims[i]=galentot;

      gahandle=NGA_Create(gadtype, ndim, dims, name2, block);

      free(name2);
      free(dims);
      free(block);

      *handle=(fortint)gahandle;
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Create should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Deallocate the array represented by handle and free any associated resources.
 *
 *  - \b GA analogous http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_destroy
 */
   void PPIDD_Destroy(fortint *handle,fortlogical *ok) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_free(mpihandle);
      else {
         int mproc=-NProcs_Work();
         mpierr=twosided_helpga_col(mproc, mpihandle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Destroy: array %d destroyed!\n",ProcID(),mpihandle);
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int ihandle = (int) *handle;
      GA_Destroy(ihandle);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Destroy should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Return the range of a distributed array held by a specified process.
 *
 *  Return <tt>.true.</tt> if successful, otherwise <tt>.false.</tt>
 *  If no array elements are owned by the process, the range is returned as [0,-1].
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_distribution
 */
   void PPIDD_Distrib(fortint *handle,fortint *rank,fortint *ilo,fortint *ihi,fortlogical *ok) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpirank=(int)*rank;
      int mpiilo;
      int mpiihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 ) {
         mpierr=mpiga_distribution( mpihandle, mpirank, &mpiilo, &mpiihi);
      }
      else {
         mpierr=twosided_helpga_distrib( mpihandle, mpirank, &mpiilo, &mpiihi);
      }

      *ilo = (fortint) mpiilo;
      *ihi = (fortint) mpiihi;
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int gahandle=(int)*handle;
      int garank=(int)*rank;
      int ndim=1;
      int *gailo;
      int *gaihi;
      int ARRAY_BASE=0;

      gailo=(int *)malloc(ndim*sizeof(int));
      gaihi=(int *)malloc(ndim*sizeof(int));
      NGA_Distribution(gahandle, garank, gailo, gaihi);
/* If no array elements are owned by process iproc, the range is returned as lo[ ]=0 and hi[ ]= -1 for all dimensions. */
      if (gailo[0]<=gaihi[0]) {
         *ilo = (fortint) (gailo[0] + 1 - ARRAY_BASE);
         *ihi = (fortint) (gaihi[0] + 1 - ARRAY_BASE);
      }
      else {
         *ilo = (fortint) (gailo[0]);
         *ihi = (fortint) (gaihi[0]);
      }
      free(gailo);
      free(gaihi);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Distrib should not be called in serial case.\n");
      exit(1);
#endif
   }



/*! Return a list of the processes that hold the data.
 *
 *  Parts of the specified patch might be actually 'owned' by several processes.
 *  np is the number of processes hold tha data (return 0  if ilo/ihi are out of bounds "0").
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_locate_region
 */
   void PPIDD_Location(fortint *handle,     /*!< array handle */
                       fortint *ilo,        /*!< lower element subscript, 1 (not 0) for fisrt element */
                       fortint *ihi,        /*!< higher element subscript */
                       fortint *map,        /*!< return start/end index for \c proclist */
                       fortint *proclist,   /*!< proc id list */
                       fortint *np,         /*!< proc number */
                       fortlogical *ok      /*!< return \c .true. if successful; otherwise \c .false. */
   ) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpisize,mpinp;
      int *mpimap,*mpiproclist;
      int mpierr,i;
      MPI_Comm mpicomm;

      mpicomm=mpigv(Compute_comm);
      MPI_Comm_size(mpicomm, &mpisize);
      mpimap=(int *)malloc(2*mpisize*sizeof(int));
      mpiproclist=(int *)malloc(mpisize*sizeof(int));

      if ( mpiga_inquire_storetype(mpihandle) == 0 ) {
         mpierr=mpiga_location( mpihandle, mpiilo, mpiihi, mpimap, mpiproclist, &mpinp);
      }
      else {
         mpierr=twosided_helpga_location( mpihandle, mpiilo, mpiihi, mpimap, mpiproclist, &mpinp);
      }

      for (i=0;i<mpinp;i++) {
         map[2*i]=(fortint)mpimap[2*i];
         map[2*i+1]=(fortint)mpimap[2*i+1];
         proclist[i]=(fortint)mpiproclist[i];
      }
      *np = (fortint) mpinp;
      free(mpimap);
      free(mpiproclist);
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int mpihandle=(int)*handle;
      int *mpiilo;
      int *mpiihi;
      int mpisize,mpinp;
      int *mpimap,*mpiproclist;
      int i;
      int ARRAY_BASE=0;
      int ndim=1;

      mpisize = GA_Nnodes();
      mpimap=(int *)malloc(2*mpisize*sizeof(int));
      mpiproclist=(int *)malloc(mpisize*sizeof(int));

      mpiilo=(int *)malloc(ndim*sizeof(int));
      mpiihi=(int *)malloc(ndim*sizeof(int));
      mpiilo[0]=(int)*ilo-1+ARRAY_BASE;
      mpiihi[0]=(int)*ihi-1+ARRAY_BASE;
      mpinp=NGA_Locate_region( mpihandle, mpiilo, mpiihi, mpimap, mpiproclist);

      for (i=0;i<mpinp;i++) {
	 map[2*i]=(fortint)(mpimap[2*i]+1-ARRAY_BASE);
	 map[2*i+1]=(fortint)(mpimap[2*i+1]+1-ARRAY_BASE);
	 proclist[i]=(fortint)mpiproclist[i];
      }
      *np = (fortint) mpinp;
      free(mpimap);
      free(mpiproclist);
      free(mpiilo);
      free(mpiihi);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Location should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Copies data from array section to the local array buffer according to starting and ending index.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_get
 */
   void PPIDD_Get(fortint *handle,fortint *ilo,fortint *ihi,void *buff,fortlogical *ok) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_get(mpihandle, mpiilo, mpiihi, buff);
      else {
         int mproc=0;
         int ielem=mpiilo;
         fortint val;
         MPI_Datatype dtype=twosided_helpga_inquire_dtype(mpihandle);
         if ( (mpiilo==mpiihi) && (dtype==MPI_INT||dtype==MPI_LONG||dtype==MPI_LONG_LONG) ) { /* PPIDD_helpga_get_inum */
            fortint nelem_valput=1;
            fortint *ibuff;

            val=twosided_helpga_one(mproc, nelem_valput, ielem, &mpihandle);
            ibuff=(fortint *)buff;
            *ibuff=(fortint)val;
         }
         else if (mpiilo <= mpiihi) { /* PPIDD_helpga_get */
            int nelem=mpiihi-mpiilo+1;

            val=twosided_helpga_extra(mproc, nelem, ielem, &mpihandle,buff);
         }
         else {
            MPIGA_Error("PPIDD_Get: starting index > ending index, handle=",mpihandle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Get: Get value from array handle= %d [%d--%d].\n",ProcID(),mpihandle,mpiilo,mpiihi);
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int mpihandle=(int)*handle;
      int ld[1]={1};
      int ndim=1;
      int *mpiilo;
      int *mpiihi;
      int ARRAY_BASE=0;

      mpiilo=(int *)malloc(ndim*sizeof(int));
      mpiihi=(int *)malloc(ndim*sizeof(int));
      mpiilo[0]=(int)*ilo-1+ARRAY_BASE;
      mpiihi[0]=(int)*ihi-1+ARRAY_BASE;
      NGA_Get(mpihandle, mpiilo, mpiihi, buff, ld);
      free(mpiilo);
      free(mpiihi);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Get should not be called in serial case.\n");
      exit(1);
#endif
   }

/*! Put local buffer data into a section of a global array according to starting and ending index.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_put
 */
   void PPIDD_Put(fortint *handle,fortint *ilo,fortint *ihi,void *buff,fortlogical *ok) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_put(mpihandle, mpiilo, mpiihi, buff);
      else {
         int mproc=-NProcs_Work();
         int ielem=mpiilo;
         MPI_Datatype dtype=twosided_helpga_inquire_dtype(mpihandle);
         if ( (mpiilo==mpiihi) && (dtype==MPI_INT||dtype==MPI_LONG||dtype==MPI_LONG_LONG) ) { /* PPIDD_helpga_put_inum */
            fortint valput;
            fortint *ibuff;

            ibuff=(fortint *)buff;
            valput=(fortint)*ibuff;
            twosided_helpga_one(mproc, valput, ielem, &mpihandle);
         }
         else if (mpiilo <= mpiihi) { /* PPIDD_helpga_put */
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra(mproc, nelem, ielem, &mpihandle,buff);
         }
         else {
            MPIGA_Error("PPIDD_Put: starting index > ending index, handle=",mpihandle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Put: Put buff numbers to array handle=%d [%d--%d].\n",ProcID(),mpihandle,mpiilo,mpiihi);
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int mpihandle=(int)*handle;
      int ld[1]={1};
      int ndim=1;
      int *mpiilo;
      int *mpiihi;
      int ARRAY_BASE=0;

      mpiilo=(int *)malloc(ndim*sizeof(int));
      mpiihi=(int *)malloc(ndim*sizeof(int));
      mpiilo[0]=(int)*ilo-1+ARRAY_BASE;
      mpiihi[0]=(int)*ihi-1+ARRAY_BASE;
      NGA_Put(mpihandle, mpiilo, mpiihi, buff, ld);
      free(mpiilo);
      free(mpiihi);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Put should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Accumulate data into a section of a global array.
 *
 * Atomic operation.  global array section (ilo, ihi) += *fac * buffer
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_acc
 */
   void PPIDD_Acc(fortint *handle,fortint *ilo,fortint *ihi,void *buff,void *fac,fortlogical *ok) {
#ifdef MPI2
      int mpihandle=(int)*handle;
      int mpiilo=(int)*ilo;
      int mpiihi=(int)*ihi;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_acc(mpihandle, mpiilo, mpiihi, buff, fac);
      else {
         if (mpiilo <= mpiihi) { /* PPIDD_helpga_acc */
            int mproc=NProcs_Work();
            int ielem=mpiilo;
            int nelem=mpiihi-mpiilo+1;

            twosided_helpga_extra_acc(mproc, nelem, ielem, &mpihandle, buff, fac);
         }
         else {
            MPIGA_Error("PPIDD_Put: starting index > ending index, handle=",mpihandle);
         }
         mpierr=0;
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Acc: Accumulate buff numbers to array handle=%d [%d--%d].\n",ProcID(),mpihandle,mpiilo,mpiihi);
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int mpihandle=(int)*handle;
      int ld[1]={1};
      int ndim=1;
      int *mpiilo;
      int *mpiihi;
      int ARRAY_BASE=0;

      mpiilo=(int *)malloc(ndim*sizeof(int));
      mpiihi=(int *)malloc(ndim*sizeof(int));
      mpiilo[0]=(int)*ilo-1+ARRAY_BASE;
      mpiihi[0]=(int)*ihi-1+ARRAY_BASE;
      NGA_Acc(mpihandle, mpiilo, mpiihi, buff, ld, fac);
      free(mpiilo);
      free(mpiihi);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Acc should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Atomically read and increment an element in an integer array.
 *
 *  Reads data from the (inum) element of a global array of integers, returns that value, and increments the (inum)
    element by a given increment. This is a fetch-and-add operation.
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_read_inc
 */
   void PPIDD_Read_inc(fortint *ihandle,fortint *inum,fortint *incr,fortint *returnval) {
#ifdef MPI2
      int mpihandle = (int) *ihandle;
      int mpiinum = (int) *inum;
      int mpiincr = (int) *incr;
      fortint mpivalue;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpivalue=mpiga_read_inc(mpihandle,mpiinum,mpiincr);
      else {                                              /* PPIDD_helpga_readinc */
         int mproc=NProcs_Work();

         mpivalue=twosided_helpga_one(mproc, *incr, mpiinum, &mpihandle);
      }
      *returnval=(fortint)mpivalue;
      if(MPIGA_Debug)printf("%5d: In PPIDD_Read_inc: fetch-and-add element[%d] of array handle=%d by increment=%d\n",
                            ProcID(),mpiinum,mpihandle,mpiincr);
#endif
#ifdef GA_TOOLS
      int handle = (int) *ihandle;
      int ndim=1;
      int *mpiinum;
      long gaincr = (long) *incr;
      long gavalue;
      int ARRAY_BASE=0;

      mpiinum=(int *)malloc(ndim*sizeof(int));
      mpiinum[0] = (int) *inum-1+ARRAY_BASE;
      gavalue=NGA_Read_inc(handle,mpiinum, gaincr);
      free(mpiinum);
      *returnval=(fortint)gavalue;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Read_inc should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Set all the elements in an array patch to zero.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_zero_patch
 */
   void PPIDD_Zero_patch(fortint *ihandle,fortint *ilo,fortint *ihi) {
#ifdef MPI2
      int mpihandle = (int) *ihandle;
      int mpiilo = (int) *ilo;
      int mpiihi = (int) *ihi;
      int mpierr=0;
      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_zero_patch(mpihandle,mpiilo,mpiihi);
      else
         MPIGA_Error("PPIDD_Zero_patch: invalid storetype, should be 0. handle=",mpihandle);

      if(mpierr!=0) MPI_Abort(mpigv(Compute_comm),911);
#endif
#ifdef GA_TOOLS
      int handle = (int) *ihandle;
      int ndim=1;
      int *mpiilo;
      int *mpiihi;
      int ARRAY_BASE=0;

      mpiilo=(int *)malloc(ndim*sizeof(int));
      mpiihi=(int *)malloc(ndim*sizeof(int));
      mpiilo[0]=(int)*ilo-1+ARRAY_BASE;
      mpiihi[0]=(int)*ihi-1+ARRAY_BASE;

      NGA_Zero_patch (handle, mpiilo, mpiihi);
      free(mpiilo);
      free(mpiihi);
#endif
   }


/*! Set all the elements of an array to zero.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_zero
 */
   void PPIDD_Zero(fortint *handle,fortlogical *ok) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      int mpierr;

      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpierr=mpiga_zero(mpihandle);
      else {
         int mproc=NProcs_Work();
         mpierr=twosided_helpga_col(mproc, mpihandle);
      }
      if(MPIGA_Debug)printf("%5d: In PPIDD_Zero: array %d has been set to zero.\n",ProcID(),mpihandle);
      if (mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int ihandle = (int) *handle;
      GA_Zero(ihandle);
      *ok = TRUE ;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Zero should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Create a new array by applying all the properties of another existing array.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_duplicate
 *  - \b MPI2  does nothing
 */
   void PPIDD_Duplicate(fortint *handlei, fortint *handlej, char *name
/*! \cond */
#if defined(FORTCL_NEXT) || defined(FORTCL_END)
	       ,fortintc lx
#endif
/*! \endcond */
      ) {
#ifdef GA_TOOLS
      int ga_a=(int)*handlei;
      int ga_b;
      char *name2, *p;
      int lxi;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
      strncpy((name2=(char *)calloc(lxi+1,1)),name,lxi);
      for (p=name2+lxi;p>=name2;p--) if (*p==' ')*p=(char)0;
      ga_b = GA_Duplicate(ga_a, name2);
      free(name2);
      *handlej=(fortint)ga_b;
#else
      printf(" ERROR: PPIDD_Duplicate should not be called in serial and MPI2 cases.\n");
      exit(1);
#endif
}


/*! Returns the name of a global array represented by the handle.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_inquire_name
 *  - \c This operation is local.
 */
   void PPIDD_Inquire_name(fortint *handle, char *name
/*! \cond */
#if defined(FORTCL_NEXT) || defined(FORTCL_END)
	       ,fortintc lx
#endif
/*! \endcond */
      ) {
#ifdef MPI2
      char *name2;
      int lxi;
      int i,len_actual;
      int mpihandle = (int) *handle;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
      if ( mpiga_inquire_storetype(mpihandle) == 0 )
         mpiga_inquire_name(mpihandle, &name2);
      else {
         twosided_helpga_inquire_name(mpihandle, &name2);
      }
      len_actual=strlen(name2);
      strncpy(name,name2,len_actual);
      for(i=len_actual;i<lxi;i++) name[i]=' ';
      if(MPIGA_Debug)printf("In PPIDD_Inquire_name: name2=%s,strlen(name2)=%d,lxi=%d\n",name2,len_actual,lxi);
#endif

#ifdef GA_TOOLS
      char *name2;
      int lxi;
      int i,len_actual;
      int gahandle = (int) *handle;

#if defined(FORTCL_NEXT) || defined(FORTCL_END)
      lxi=(int)lx;
#ifdef FORTINTC_DIVIDE
      lxi=(int)lx/FORTINTC_DIVIDE;
#endif
#else
      lxi=strlen(name);
#endif
/*      strcpy(name2=malloc(80*sizeof(char)),GA_Inquire_name(gahandle));
      strncpy(name,name2,strlen(name2)); */
      name2=GA_Inquire_name(gahandle);
      len_actual=strlen(name2);
      strncpy(name,name2,len_actual);
      for(i=len_actual;i<lxi;i++) name[i]=' ';
      if(MPIGA_Debug)printf("In PPIDD_Inquire_name: name2=%s,strlen(name2)=%d,lxi=%d\n",name2,len_actual,lxi);
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Inquire_name should not be called in serial case.\n");
      exit(1);
#endif
}


/*! Returns the storetype of a global array represented by the handle.
 *
 *     storetype: number of helper processes for storing this global array.
 *  - \c storetype=0 : Normal distributed array stored across the distributed processes
 *  - \c storetype>=1: Low-latency array stored on one or more helper processes (effective only when helper process is enabled).
 *  - \c This operation is local.
 */
   void PPIDD_Inquire_stype(fortint *handle, fortint *storetype) {
#ifdef MPI2
      int mpihandle = (int) *handle;
      *storetype=(fortint)mpiga_inquire_storetype(mpihandle);
#endif
#ifdef GA_TOOLS
      *storetype=(fortint)0;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      printf(" ERROR: PPIDD_Inquire_stype should not be called in serial case.\n");
      exit(1);
#endif
   }


/*! Get the amount of memory (in bytes) used in the allocated distributed arrays on the calling processor.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_inquire_memory
 */
   void PPIDD_Inquire_mem(fortint *mem_used) {
#ifdef MPI2
      long localmem;
      localmem=mpiga_localmem();
      *mem_used=(fortint)localmem;
#endif
#ifdef GA_TOOLS
      size_t localmem;
      localmem=GA_Inquire_memory();
      *mem_used=(fortint)localmem;
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *mem_used=(fortint)0;
#endif
   }



/*! Create a set of mutex variables that can be used for mutual exclusion.
 *
 *  Returns <tt>.true.</tt> if the operation succeeded or <tt>.false.</tt> when failed. It is a collective operation.

 *  - For \b GA, \c storetype doesn't take effect, analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_create_mutexes
 *  - For \b MPI2, there are two kinds of mutexes:
 *   - (1) \c storetype=0: mutex data will be stored by a global array across the distributed processes.
 *        Due to the poor performance of MPI-2 one-sided communication, this kind of mutex is very slow.
 *   - (2) \c storetype=1: mutex data will be stored on the helper process. Improve the performance of (1), and all communications are based on two-sided opeartions.
 *
 *   When helper process is disabled, \c storetype doesn't take effect, and mutex data are always stored across the distributed processes.

 */
   void PPIDD_Create_mutexes(fortint *storetype,fortint *number,fortlogical *ok) {
#ifdef MPI2
      int stype     = (int) *storetype;
      int mpinumber = (int) *number;
      int mpierr;

      if (use_helper_server==0) {
        mpierr=mpiga_create_mutexes(mpinumber);      /* mutexes data store by a global array across the distributed processes */
      }
      else {
        if (stype==0)
         mpierr=mpiga_create_mutexes(mpinumber);      /* mutexes data store by a global array across the distributed processes */
        else
         mpierr=alloc_general_helpmutexes(mpinumber);  /* mutexes data on helper process */
      }

      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int mpinumber = (int) *number;
      int mpierr;
      mpierr=GA_Create_mutexes(mpinumber);
/* This is one of exceptions in GA (see global/src/capi.c) : Returns [1] if the operation succeeded or [0] when failed */
      if(mpierr==1) *ok = TRUE ;
      else *ok = FALSE ;
      if(MPIGA_Debug)printf("In PPIDD_Create_Mutexes: mpierr=%d, ok=%d.\n",mpierr,(int)*ok);
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *ok = TRUE ;
#endif
   }


/*! Lock a mutex object identified by a given mutex number.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_lock
 */
   void PPIDD_Lock_mutex(fortint *inum) {
#ifdef MPI2
      int mpiinum = (int) *inum;
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_lock_mutex(mpiinum);   /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=lock_general_helpmutex(mpiinum); /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpigv(Compute_comm),911);
#endif
#ifdef GA_TOOLS
      int mpiinum = (int) *inum;
      GA_Lock(mpiinum);
#endif
   }


/*! Unlock  a mutex object identified by a given mutex number.
 *
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_unlock
 */
   void PPIDD_Unlock_mutex(fortint *inum) {
#ifdef MPI2
      int mpiinum = (int) *inum;
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_unlock_mutex(mpiinum);    /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=unlock_general_helpmutex(mpiinum);  /* mutexes data on helper process */
      if(mpierr!=0) MPI_Abort(mpigv(Compute_comm),911);
#endif
#ifdef GA_TOOLS
      int mpiinum = (int) *inum;
      GA_Unlock(mpiinum);
#endif
   }


/*! Destroy the set of mutexes created with PPIDD_Create_mutexes.
 *
 * Returns <tt>.true.</tt> if the operation succeeded or <tt>.false.</tt> when failed. This is a collective operation.
 *  - \b GA analogous to http://www.emsl.pnl.gov/docs/global/ga_ops.html#ga_destroy_mutexes
 */
   void PPIDD_Destroy_mutexes(fortlogical *ok) {
#ifdef MPI2
      int mpierr;
      if ( mpigv(nmutex) > 0 )
         mpierr=mpiga_destroy_mutexes();  /* mutexes data store by a global array across the distributed processes */
      else
         mpierr=free_general_helpmutexes();   /* mutexes data on helper process */
      if(mpierr==0) *ok = TRUE ;
      else *ok = FALSE ;
#endif
#ifdef GA_TOOLS
      int mpierr;
      mpierr=GA_Destroy_mutexes();
/* This is one of exceptions in GA (see global/src/capi.c) : Returns [1] if the operation succeeded or [0] when failed */
      if(mpierr==1) *ok = TRUE ;
      else *ok = FALSE ;
      if(MPIGA_Debug)printf("In PPIDD_Destroy_Mutexes: mpierr=%d, ok=%d.\n",mpierr,(int)*ok);
#endif
#if !defined(MPI2) && !defined(GA_TOOLS)
      *ok = TRUE ;
#endif
   }
