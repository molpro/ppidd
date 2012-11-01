
/*! \mainpage PPIDD Reference Manual
 *
 * Parallel Programming Interface for Distributed Data (PPIDD) Reference
 * Manual. The Fortran interface subroutine descriptions can be found
 * in ppidd_fortran.c and ppidd_share.h, and the C interface subroutine
 * descriptions can be found in ppidd_c.c and ppidd_share.h.

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
       |   |-- mpi_helpmutex.c                Mutex source file using helper process
       |   |-- mpi_nxtval.c                   NXTVAL source file
       |   |-- mpi_nxtval.h                   NXTVAL header file
       |   |-- mpi_utils.c                    MPI utility source file
       |   |-- mpi_utils.h                    MPI utility header file
       |   |-- mpiga_base.c                   Source code for distributed data structure
       |   |-- mpiga_base.h                   Head file for distributed data structure
       |   |-- mpimutex-hybrid.c              Mutex source file using distributed processes
       |   |-- mpimutex.h                     Mutex header file using distributed processes
       |   |-- ppidd_c.c                      C interface source code
       |   |-- ppidd_c.h                      C interface header file
       |   |-- ppidd_doxygen.h                PPIDD document main page file
       |   |-- ppidd_dtype.h                  PPIDD data type header file
       |   |-- ppidd_eaf_c.c                  C interface source code for EAF
       |   |-- ppidd_eaf_c.h                  C interface header file for EAF
       |   |-- ppidd_eaf_fortran.c            Fortran interface source code for EAF
       |   |-- ppidd_eaf_fortran.h            Fortran interface header file for EAF
       |   |-- ppidd_eaf_share.h              Code shared by both Fortran and C interfaces for EAF
       |   |-- ppidd_fortran.c                Fortran interface source code for PPIDD
       |   |-- ppidd_fortran.h                Fortran interface header file for PPIDD
       |   |-- ppidd_sf_c.c                   C interface source code for SF
       |   |-- ppidd_sf_c.h                   C interface header file for SF
       |   |-- ppidd_sf_fortran.c             Fortran interface source code for SF
       |   |-- ppidd_sf_fortran.h             Fortran interface header file for SF
       |   |-- ppidd_sf_share.h               Code shared by both Fortran and C interfaces for SF
       |   |-- ppidd_share.h                  Code shared by both Fortran and C interfaces
       |   `-- ppidd_undefdtype.h             Fortran data type header file for C interface
       `-- ./test                         Test code directory for PPIDD library
           |-- GNUmakefile                    Makefile for test directory
           |-- ppidd_ctest.c                  C test program
           |-- ppidd_test.F                   Fortran test program
           |-- sizeofctypes.c                 Code for determining the size of C data types
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
     make MPICC=mpiicc MPIFC=mpiifort _I8_=y FFLAGS='-i8'

     or
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8' INCLUDE=/software/intel/mpi/4.0.0.025/intel64/include \
     MPILIB='-L/software/intel/mpi/4.0.0.025/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /software/intel/mpi/4.0.0.025/intel64/lib -Xlinker -rpath -Xlinker \
     /opt/intel/mpi-rt/4.0.0 -lmpi -lmpigf -lmpigi -lpthread -lpthread -lpthread -lpthread -lrt -L/usr/lib64 -libverbs -lm'
     (the front part for MPILIB option comes from `mpiifort -show`, and the rear part '-L/usr/lib64 -libverbs -lm' is used to link with Infiniband network.)

     2. Build Global Arrays version of PPIDD.
     Global Arrays should be installed prior to building PPIDD. As mentioned in Global Arrays documentation(http://www.emsl.pnl.gov/docs/global),
     there are three possible ways for building GA: (1) GA with MPI; (2) GA with TCGMSG-MPI; and (3) GA with TCGMSG.
     PPIDD can be built with either of these interfaces. As the structure of GA library has changed significantly since version 5.0,
     the building methods of PPIDD differ correspondingly.

     A. For GA version lower than 5.0:
     (1) GA with MPI:
     make MPICC=mpiicc MPIFC=mpiifort _I8_=y FFLAGS='-i8 -Vaxlib' INCLUDE='../../../ga-4-3-3/include /software/intel/mpi/4.0.0.025/intel64/include' MPILIB='-L/usr/lib64 -libverbs -lm'

     or
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8 -Vaxlib' INCLUDE='../../../ga-4-3-3/include /software/intel/mpi/4.0.0.025/intel64/include' \
     MPILIB='-L/software/intel/mpi/4.0.0.025/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /software/intel/mpi/4.0.0.025/intel64/lib -Xlinker -rpath -Xlinker \
     /opt/intel/mpi-rt/4.0.0 -lmpi -lmpigf -lmpigi -lpthread -lpthread -lpthread -lpthread -lrt -L/usr/lib64 -libverbs -lm'


     (2) GA with TCGMSG-MPI:
     make MPICC=mpiicc MPIFC=mpiifort _I8_=y FFLAGS='-i8 -Vaxlib' INCLUDE=../../../ga-4-3-3/include MPILIB='-L/usr/lib64 -libverbs -lm'

     or
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8 -Vaxlib' INCLUDE=../../../ga-4-3-3/include \
     MPILIB='-L/software/intel/mpi/4.0.0.025/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /software/intel/mpi/4.0.0.025/intel64/lib -Xlinker -rpath -Xlinker \
     /opt/intel/mpi-rt/4.0.0 -lmpi -lmpigf -lmpigi -lpthread -lpthread -lpthread -lpthread -lrt -L/usr/lib64 -libverbs -lm'

     (3) GA with TCGMSG:
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8 -Vaxlib' INCLUDE=../../../ga-4-3-3/include MPILIB=../../../ga-4-3-3/lib/LINUX64/lib/

     B. For GA version higher than (or equal to) 5.0, BUILD=GA_TCGMSG, GA_TCGMSGMPI or GA_MPI must be specified explicitly:
     (1) GA with MPI:
     make MPICC=mpiicc MPIFC=mpiifort _I8_=y FFLAGS='-i8 -Vaxlib' BUILD=GA_MPI INCLUDE='../../../ga-5-0-2-install/include /software/intel/mpi/4.0.0.025/intel64/include' \
     MPILIB='-L/usr/lib64 -libverbs -lm -L/software/intel/mkl/10.0.1.014/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core'

     or
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8 -Vaxlib' BUILD=GA_MPI INCLUDE='../../../ga-5-0-2-install/include /software/intel/mpi/4.0.0.025/intel64/include' \
     MPILIB='-L/software/intel/mpi/4.0.0.025/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /software/intel/mpi/4.0.0.025/intel64/lib -Xlinker -rpath -Xlinker \
     /opt/intel/mpi-rt/4.0.0 -lmpi -lmpigf -lmpigi -lpthread -lpthread -lpthread -lpthread -lrt -L/usr/lib64 -libverbs -lm \
     -L/software/intel/mkl/10.0.1.014/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core'
     (the front part for MPILIB option comes from `mpiifort -show`, the middle part '-L/usr/lib64 -libverbs -lm' is used to link with Infiniband network, and the real part
     -L/software/intel/mkl/10.0.1.014/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core is used to link with blas library if needed.)

     (2) GA with TCGMSG-MPI:
     make MPICC=mpiicc MPIFC=mpiifort _I8_=y FFLAGS='-i8 -Vaxlib' BUILD=GA_TCGMSGMPI INCLUDE=../../../ga-5-0-2-install/include \
     MPILIB='-L/usr/lib64 -libverbs -lm -L/software/intel/mkl/10.0.1.014/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core'

     or
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8 -Vaxlib' BUILD=GA_TCGMSGMPI INCLUDE=../../../ga-5-0-2-install/include \
     MPILIB='-L/software/intel/mpi/4.0.0.025/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /software/intel/mpi/4.0.0.025/intel64/lib -Xlinker -rpath -Xlinker \
     /opt/intel/mpi-rt/4.0.0 -lmpi -lmpigf -lmpigi -lpthread -lpthread -lpthread -lpthread -lrt -L/usr/lib64 -libverbs -lm \
     -L/software/intel/mkl/10.0.1.014/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core'

     (3) GA with TCGMSG:
     make CC=icc FC=ifort _I8_=y FFLAGS='-i8 -Vaxlib' BUILD=GA_TCGMSG INCLUDE=../../../ga-5-0-2-install/include \
     MPILIB='-L/software/intel/mkl/10.0.1.014/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core'

  </pre>
 *
 */
