
/*! \mainpage PPIDD Reference Manual
 *
 * Parallel Programming Interface for Distributed Data (PPIDD) Reference
 * Manual. The Fortran interface subroutine descriptions can be found
 * in ppidd_fortran.cpp and ppidd_share.h, and the C interface subroutine
 * descriptions can be found in ppidd_c.cpp and ppidd_share.h.

//   Directory and file structure in PPIDD:
  <pre>

       .                                  PPIDD root directory
       |-- GNUmakefile                    Makefile to build PPIDD
       |-- README                         README for PPIDD
       |-- ./doc                          Document directory
       |   `-- Doxyfile                       Document configuration file
       |-- ./lib                          Final location of PPIDD library
       |-- ./src                          Source code directory for PPIDD library
       |   |-- GNUmakefile                    Makefile for src directory
       |   |-- ppidd_machines.h               Head file for machine-related settings
       |   |-- mpi_helpmutex.cpp              Mutex source file using helper process
       |   |-- mpi_nxtval.cpp                 NXTVAL source file
       |   |-- mpi_nxtval.h                   NXTVAL header file
       |   |-- mpi_utils.cpp                  MPI utility source file
       |   |-- mpi_utils.h                    MPI utility header file
       |   |-- mpiga_base.cpp                 Source code for distributed data structure
       |   |-- mpiga_base.h                   Head file for distributed data structure
       |   |-- mpimutex-hybrid.cpp            Mutex source file using distributed processes
       |   |-- mpimutex.h                     Mutex header file using distributed processes
       |   |-- ppidd_c.cpp                    C interface source code
       |   |-- ppidd_c.h                      C interface header file
       |   |-- ppidd_doxygen.h                PPIDD document main page file
       |   |-- ppidd_dtype.h                  PPIDD data type header file
       |   |-- ppidd_eaf_c.cpp                C interface source code for EAF
       |   |-- ppidd_eaf_c.h                  C interface header file for EAF
       |   |-- ppidd_eaf_fortran.cpp          Fortran interface source code for EAF
       |   |-- ppidd_eaf_fortran.h            Fortran interface header file for EAF
       |   |-- ppidd_eaf_share.h              Code shared by both Fortran and C interfaces for EAF
       |   |-- ppidd_fortran.cpp              Fortran interface source code for PPIDD
       |   |-- ppidd_fortran.h                Fortran interface header file for PPIDD
       |   |-- ppidd_sf_c.cpp                 C interface source code for SF
       |   |-- ppidd_sf_c.h                   C interface header file for SF
       |   |-- ppidd_sf_fortran.cpp           Fortran interface source code for SF
       |   |-- ppidd_sf_fortran.h             Fortran interface header file for SF
       |   |-- ppidd_sf_share.h               Code shared by both Fortran and C interfaces for SF
       |   |-- ppidd_share.h                  Code shared by both Fortran and C interfaces
       |   `-- ppidd_undefdtype.h             Fortran data type header file for C interface
       `-- ./test                         Test code directory for PPIDD library
           |-- GNUmakefile                    Makefile for test directory
           |-- ppidd_ctest.cpp                C test program
           |-- ppidd_test.F                   Fortran test program
           |-- sizeofctypes.cpp               Code for determining the size of C data types
           |-- sizeoffortypes.F               Code for determining the size of Fortran data types
           |-- ppidd_test.F                   Fortran test program
           `-- ppidd_test.out                 Output of test example

  </pre>

//   Some examples of building the PPIDD library:
 *    The following examples are tested on a x86_64//Linux machine on which Intel Fortran and C compilers,
 *    MPI2-aware Intel Fortran and C compilers, and Intel MPI library are available.
 *    Please be aware the options might be different on other machines.
  <pre>

     1. Build MPI-2 version of PPIDD:
     ./configure --with-mpi2

     2. Build Global Arrays version of PPIDD.
     ./configure --with-ga=/path/to/ga-install/include

  </pre>
 *
 */