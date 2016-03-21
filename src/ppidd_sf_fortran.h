
/* ------------------------------------------------------------- *\
   FORTRAN Interface of PPIDD Shared Files Library
   ( External PPIDD SF routines which are used to be called by FORTRAN routines)
   =====================================================================

   All PPIDD subroutines are written in C.  PPIDD subroutines called
   from FORTRAN go through a C wrapper.  All FORTRAN wrapper
   subroutines are prefixed with PPIDD_.

\* ------------------------------------------------------------- */

#ifndef __PPIDD_SF_FORTRAN_H__
#define __PPIDD_SF_FORTRAN_H__

 #     define PPIDD_Sf_create               FC_FUNC_(ppidd_sf_create,PPIDD_SF_CREATE)
 #     define PPIDD_Sf_write                FC_FUNC_(ppidd_sf_write,PPIDD_SF_WRITE)
 #     define PPIDD_Sf_read                 FC_FUNC_(ppidd_sf_read,PPIDD_SF_READ)
 #     define PPIDD_Sf_wait                 FC_FUNC_(ppidd_sf_wait,PPIDD_SF_WAIT)
 #     define PPIDD_Sf_waitall              FC_FUNC_(ppidd_sf_waitall,PPIDD_SF_WAITALL)
 #     define PPIDD_Sf_destroy              FC_FUNC_(ppidd_sf_destroy,PPIDD_SF_DESTROY)
 #     define PPIDD_Sf_errmsg               FC_FUNC_(ppidd_sf_errmsg,PPIDD_SF_ERRMSG)

#endif
