
/*! \file
  Overwrite data types (fortint and forlogical) for PPIDD C interface.
*/

#ifndef __PPIDD_UNDEFDTYPE_H__
#define __PPIDD_UNDEFDTYPE_H__

#ifdef fortint
#undef fortint
#endif
#define fortint int64_t

#ifdef fortlogical
#undef fortlogical
#endif
#define fortlogical int64_t

#ifdef FORTCL_NEXT
#undef FORTCL_NEXT
#endif

#ifdef FORTCL_END
#undef FORTCL_END
#endif

#endif
