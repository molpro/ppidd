
/* ------------------------------------------------------------- *\
   FORTRAN Interface of PPIDD Exclusive Access File Library
   ( External PPIDD EAF routines which are used to be called by FORTRAN routines)
   =====================================================================

   All PPIDD subroutines are written in C.  PPIDD subroutines called
   from FORTRAN go through a C wrapper.  All FORTRAN wrapper
   subroutines are prefixed with PPIDD_.

   The following C definitions substitute the FORTRAN wrapper
   subroutine name with the correct machine dependent FORTRAN
   external name.

   Note: FORTRAN externals are generally all lowercase, but may
   be uppercase.  See ppidd_machines.h for details.

\* ------------------------------------------------------------- */

#ifndef __PPIDD_EAF_FORTRAN_H__
#define __PPIDD_EAF_FORTRAN_H__

 # include "ppidd_machines.h"

 #     define PPIDD_Eaf_open                 FC_FUNC_(ppidd_eaf_open,PPIDD_EAF_OPEN)
 #     define PPIDD_Eaf_write                FC_FUNC_(ppidd_eaf_write,PPIDD_EAF_WRITE)
 #     define PPIDD_Eaf_awrite               FC_FUNC_(ppidd_eaf_awrite,PPIDD_EAF_AWRITE)
 #     define PPIDD_Eaf_read                 FC_FUNC_(ppidd_eaf_read,PPIDD_EAF_READ)
 #     define PPIDD_Eaf_aread                FC_FUNC_(ppidd_eaf_aread,PPIDD_EAF_AREAD)
 #     define PPIDD_Eaf_wait                 FC_FUNC_(ppidd_eaf_wait,PPIDD_EAF_WAIT)
 #     define PPIDD_Eaf_waitall              FC_FUNC_(ppidd_eaf_waitall,PPIDD_EAF_WAITALL)
 #     define PPIDD_Eaf_probe                FC_FUNC_(ppidd_eaf_probe,PPIDD_EAF_PROBE)
 #     define PPIDD_Eaf_close                FC_FUNC_(ppidd_eaf_close,PPIDD_EAF_CLOSE)
 #     define PPIDD_Eaf_delete               FC_FUNC_(ppidd_eaf_delete,PPIDD_EAF_DELETE)
 #     define PPIDD_Eaf_length               FC_FUNC_(ppidd_eaf_length,PPIDD_EAF_LENGTH)
 #     define PPIDD_Eaf_truncate             FC_FUNC_(ppidd_eaf_truncate,PPIDD_EAF_TRUNCATE)
 #     define PPIDD_Eaf_errmsg               FC_FUNC_(ppidd_eaf_errmsg,PPIDD_EAF_ERRMSG)

#endif
