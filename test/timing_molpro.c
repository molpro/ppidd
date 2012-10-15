#include "util/machines.h"

#define	SECOND        FORT_Extern(second,SECOND)
#define	TIMING_MOLPRO FORT_Extern(timing_molpro,TIMING_MOLPRO)
#define WALLCL        FORT_Extern(wallcl,WALLCL)

#ifdef _WIN32
int gettimeofday (struct timeval *tp, int *tzp) {
return 0;
}
#endif

#ifdef __cplusplus
extern "C" {void TIMING_MOLPRO(double&, double&, double&);}
#endif

#ifdef _WIN32
static clock_t itt;
#else
static clock_t it; struct tms itt;
#endif

double SECOND()
{
#if defined(CATAMOUNT) || defined(_WIN32)
 return (double)0;
#else
 it=times(&itt);
 return (double) itt.tms_utime / (double) TIMER_TICK;
#endif
}

double WALLCL()
{
#if defined(sunx) || defined(CATAMOUNT) || defined(_WIN32) /* actually this should probably work everywhere */
 return (double) time(NULL);
#else
 it=times(&itt);
 return (double) it / (double) TIMER_TICK;
#endif
}

void TIMING_MOLPRO(double *cpu, double *sys, double *wall)
{
#if defined(CATAMOUNT) || defined(_WIN32)
 *wall=*cpu=*sys= (double)0;
#else
 it=times(&itt);
 *cpu=(double) itt.tms_utime / (double) TIMER_TICK; /* same as SECOND() */
 *sys=(double) itt.tms_stime / (double) TIMER_TICK;
 *wall=(double) it / (double) TIMER_TICK; /* same as WALLCL() */
#endif
}
