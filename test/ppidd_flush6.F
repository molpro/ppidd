      subroutine ppidd_flush6
! Use FLUSH statement instead of the FLUSH intrinsic for the Fortran 2003 standard
#ifdef PPIDD_f2003
      flush(6)
#endif
! One can use the FLUSH intrinsic routine if the compiler support it.
! Some compilers(eg. NAG) might not provide such intrinsic.
! When F2003 Fortran standard is fully supported for the entire PPIDD,
! the following line can be removed.
      call flush(6) 
      end
