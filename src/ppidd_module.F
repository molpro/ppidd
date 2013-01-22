      module ppidd

      interface
        subroutine ppidd_initialize(nargs, args)
        integer, intent(in) :: nargs
        character(len=*),dimension(nargs),intent(in) :: args
        end subroutine ppidd_initialize
      end interface

      end module ppidd

      subroutine ppidd_initialize(nargs, args)
        use, intrinsic :: iso_c_binding, only: c_ptr,c_loc,c_null_char,
#ifdef _I8_
     >                                         c_long_long
#else
     >                                         c_int
#endif
        implicit none
        interface
          subroutine ppidd_initialize_f2c(nargs, args)
     >               bind(c,name="ppidd_initialize_f2c")
            use, intrinsic :: iso_c_binding, only: c_ptr,
#ifdef _I8_
     >                                         c_long_long
#else
     >                                         c_int
#endif
            implicit none
#ifdef _I8_
            integer(c_long_long), value :: nargs
#else
            integer(c_int), value :: nargs
#endif
            type(c_ptr), dimension(nargs) :: args
          end subroutine ppidd_initialize_f2c
        end interface
        integer, intent(in) :: nargs
        character(len=*),dimension(nargs),intent(in) :: args
        integer, parameter  :: lenmax=1024
        character(len=lenmax),dimension(nargs),target :: array_strings
#ifdef _I8_
        integer(c_long_long) :: nstrs
#else
        integer(c_int) :: nstrs
#endif
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
