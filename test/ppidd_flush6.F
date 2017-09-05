! Use FLUSH statement instead of the FLUSH intrinsic.
! Since the Fortran 2003 standard, there is a FLUSH statement
! that should be preferred over the FLUSH intrinsic.
      subroutine ppidd_flush6
      flush(6)
      end
