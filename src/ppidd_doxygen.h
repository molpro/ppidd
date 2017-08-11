
/*! \mainpage PPIDD Reference Manual
 *
 * Parallel Programming Interface for Distributed Data (PPIDD) Reference
 * Manual. The function descriptions can be found in ppidd.cpp.

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
       |   |-- mpi_helpmutex.cpp              Mutex source file using helper process
       |   |-- mpi_nxtval.cpp                 NXTVAL source file
       |   |-- mpi_nxtval.h                   NXTVAL header file
       |   |-- mpi_utils.cpp                  MPI utility source file
       |   |-- mpi_utils.h                    MPI utility header file
       |   |-- mpiga_base.cpp                 Source code for distributed data structure
       |   |-- mpiga_base.h                   Head file for distributed data structure
       |   |-- mpimutex-hybrid.cpp            Mutex source file using distributed processes
       |   |-- mpimutex.h                     Mutex header file using distributed processes
       |   |-- ppidd.cpp                      PPIDD interface source code
       |   |-- ppidd.h                        C interface header file
       |   |-- ppidd_doxygen.h                PPIDD document main page file
       |   |-- ppidd_eaf.cpp                  PPIDD interface source code for EAF
       |   |-- ppidd_eaf.h                    C interface header file for EAF
       |   |-- ppidd_eaf_share.h              Code shared by both Fortran and C interfaces for EAF
       |   |-- ppidd_sf.cpp                   PPIDD interface source code for SF
       |   |-- ppidd_sf.h                     C interface header file for SF
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
