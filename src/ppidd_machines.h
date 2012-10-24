
/*! \file
 * Include file for all PPIDD C files
 *
 * This include file should always be included by all PPIDD C files.
 *
 * \b FORT_Extern(subroutine_name,SUBROUTINE_NAME)
 *
 * Description: C Macro that converts the name of a C subroutine
 * into an external name that can be called by FORTRAN 77/90/95.
 *
 * This is commonly used for subroutines that need to be called
 * from both FORTRAN and C.  The subroutine is written in C and
 * a FORTRAN wrapper function is created to call the C code.
 *
 * FORTRAN externals are machine dependent.  Subroutine objects
 * names are generally augmented with 0, 1 or 2 underscores.
 *
 * _UNDERSCORES should be defined with the number of underscores
 * required by FORTRAN on a particular machine. If this is unknown,
 * compile a simple FORTRAN subroutine with the -c option and
 * use 'nm' to look at the object file.
 *
 * This header file is a subset of machines.h, and mainly includes:
 * (1) system include files on different machines
 * (2) definitions of C data types in some special cases
 * (3) definitions of LARGEFILE settings for Parallel IO (in GA)
 * (4) definitions of FORTCL settings for passing character string between Fortran and C
 * (5) definitions of Fortran data types in C
 * (6) definitions of _UNDERSCORES and Macro FORT_Extern for objects which can be called in Fortran
 */

/* PPIDD header file for C programs:
 * include files
 * default definition values, override with -Dname=[definition] */

/* cpp flag
 * _AIX		AIX
 * __APPLE__	Darwin
 * __CYGWIN__	Cygwin
 * __hpux	HP-UX
 * sgi		IRIX64
 * linux	Linux
 * sun		SunOS
 * SX		SUPER-UX
 * __uxp__	Fujitsu
 * _WIN32	Windows
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
#define ssize_t int
#define mode_t unsigned long
#include <io.h>
#define F_OK 0
#else
#include <unistd.h>
#ifndef _AIX
#include <sys/unistd.h>
#endif
#endif

#ifdef SX
#include <stdint.h>
#include <sys/socket.h> /* for gethostname */
#endif


/* definitions of C data types in some special cases */
#if ( defined(_WIN32) && !defined(__cplusplus) ) || defined(__uxp__)
typedef char int8_t;
#endif

#if defined(_WIN32) || defined(__uxp__)
typedef short int16_t;
typedef int int32_t;
#endif


/* definitions of LARGEFILE settings for Parallel IO (in GA) */

#if defined(linux) && !defined(NOLARGEFILES)
#if defined(__x86_64__) || defined(__ia64)
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
/* #define _USE_LARGEFILE64  Probably this is needed too, but it seems to work!? */
#endif
#else
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#define _FILE_OFFSET_BITS 64
#endif
#endif

/* definitions of FORTCL settings for passing character string between Fortran and C */

#if defined(_AIX) || defined(__APPLE__) || defined( __CYGWIN__) || defined(__hpux) || defined(linux) || defined(sgi) || defined(sun) || defined(SX) || defined(__uxp__) || defined(_WIN32)
#define FORTCL_END
#endif

#if defined(FORTCL_END) || defined(FORTCL_NEXT)
#define FORTCL
#endif


/* definitions of Fortran data types in C */

#ifndef FORTINT
#ifdef _I8_
#if defined(_AIX) || defined( __CYGWIN__) || defined(__hpux) || defined(linux) || defined(sgi)
#define FORTINT long long
/* set macro EXT_INT64 and EXT_INT in order to keep consistent with GA */
#define EXT_INT64
#define EXT_INT
#endif
#else
#if defined(sgi) || defined(__APPLE__) || defined(SX) || defined(__uxp__) || ( defined(__x86_64__) && defined(linux) )
#define FORTINT int
#endif
#endif
#endif
#ifndef FORTINT
#define FORTINT	long
#define EXT_INT
#endif
typedef FORTINT fortint ; /* fortran integer type */

#ifndef FORTINTC
#ifdef __uxp__
#define FORTINTC int
#endif
#ifdef _I8_
#if defined(_AIX) || defined( __CYGWIN__) || defined(__hpux) || defined(linux) || defined(sgi)
#define FORTINTC long
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


/* definitions of _UNDERSCORES and FORT_Extern for objects which can be called in Fortran */

#if defined(NAME_LU) || defined(_AIX) || defined(__APPLE__) || defined( __CYGWIN__) || defined(__hpux) || defined(linux) || defined(sgi) || defined(sun) || defined(SX) || defined(_WIN32)
#define _UNDERSCORES 1
#endif
#if defined(NAME_L)
#define _UNDERSCORES 0
#endif

#if defined _UNDERSCORES

#ifdef FORT_UPPERCASE

#if _UNDERSCORES == 0
#define FORT_Extern(funct,FUNCT) FUNCT
#elif _UNDERSCORES == 1
#define FORT_Extern(funct,FUNCT) FUNCT ## _
#elif _UNDERSCORES == 2
#define FORT_Extern(funct,FUNCT) FUNCT ## __
#else
#error "_UNDERSCORES not properly defined."
#endif

#else

#if _UNDERSCORES == 0
#define FORT_Extern(funct,FUNCT) funct
#elif _UNDERSCORES == 1
#define FORT_Extern(funct,FUNCT) funct ## _
#elif _UNDERSCORES == 2
#define FORT_Extern(funct,FUNCT) funct ## __
#else
#error "_UNDERSCORES not properly defined."
#endif

#endif

#else
#error "_UNDERSCORES not defined."
#endif

#endif /*  __PPIDD_MACHINES_H__ */
