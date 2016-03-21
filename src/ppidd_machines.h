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

#ifndef NOLARGEFILES
#define _FILE_OFFSET_BITS 64
#endif

#endif /*  __PPIDD_MACHINES_H__ */
