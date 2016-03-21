/*! \file
 * Include file for all PPIDD C files
 *
 * This include file should always be included by all PPIDD C files.
 *
 * This header file is a subset of machines.h, and mainly includes:
 * (1) system include files on different machines
 * (2) definitions of C data types in some special cases
 * (3) definitions of LARGEFILE settings for Parallel IO (in GA)
 * (4) definitions of FORTCL settings for passing character string between Fortran and C
 * (5) definitions of Fortran data types in C
 */

#ifndef __PPIDD_MACHINES_H__
#define __PPIDD_MACHINES_H__

#if defined(__linux) && !defined(linux) /* needed for Linux/ppc64 */
#define linux
#endif

/* system include files on different machines */

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef _WIN32
#include <sys/param.h>
#endif
#include <sys/stat.h>
#include <time.h>
#ifndef _WIN32
#include <sys/time.h>
#include <sys/times.h>
#endif
#include <sys/types.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#ifndef _AIX
#include <sys/unistd.h>
#endif
#endif

#ifdef SX
#include <stdint.h>
#endif

/* definitions of LARGEFILE settings for Parallel IO (in GA) */

#ifndef NOLARGEFILES
#define _FILE_OFFSET_BITS 64
#endif

/* every current supported system is FORTCL_END, some historical ones are FORTCL_NEXT */
#define FORTCL_END

/* definitions of Fortran data types in C */

#ifndef FORTINTC
#ifdef _I8_
#if defined(_AIX) || defined(__hpux) || defined(linux) || defined(sgi)
#define FORTINTC long
#elif defined(__CYGWIN__)
#define FORTINTC int
#endif
#else
#if defined(sgi)
#define FORTINTC long
#endif
#endif
#endif
#ifndef FORTINTC
#define FORTINTC FORTINT
#endif
typedef FORTINTC fortintc ; /* fortran character string length type */

#endif /*  __PPIDD_MACHINES_H__ */
