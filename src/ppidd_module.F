      module ppidd

      interface
        subroutine ppidd_initialize_f2c(nargs, args)
     >             bind(c,name="ppidd_initialize_f2c")
          use, intrinsic :: iso_c_binding, only: c_int, c_ptr
          implicit none
          integer(c_int), value :: nargs
          type(c_ptr), dimension(nargs) :: args
        end subroutine ppidd_initialize_f2c
      end interface

      contains

      subroutine ppidd_initialize(nargs, args)
        use, intrinsic :: iso_c_binding, only: c_int,c_ptr,c_loc,
     >                    c_null_char
        implicit none
        integer, intent(in) :: nargs
        character(len=*),dimension(nargs),intent(in) :: args
        integer, parameter  :: lenmax=1024
        character(len=lenmax),dimension(nargs),target :: array_strings
        integer(c_int) :: nstrs
        type(c_ptr), dimension(nargs) :: pt_string
        integer i,l

        nstrs=nargs
        do i = 1, nargs
           l=len_trim(args(i))
           if (l.eq.0) stop ! need to call error routine, argument length is zero, likely nargs is wrong
           if (l.ge.lenmax) stop ! need to call error routine, argument too long
           array_strings(i) = args(i)(1:l)//c_null_char
           pt_string(i) = c_loc(array_strings(i))
        end do

        call ppidd_initialize_f2c(nstrs,pt_string)

      end subroutine ppidd_initialize
  
      end module ppidd
