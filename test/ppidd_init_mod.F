      module ppidd_init_mod

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
     >                    c_null_ptr,c_null_char
        implicit none
        integer, intent(in) :: nargs
        integer, parameter  :: lenmax=1024
        character(len=lenmax),dimension(nargs),intent(in) :: args
        character(len=lenmax),dimension(nargs),target :: array_strings
        integer(c_int) :: nstrs
        type(c_ptr), dimension(nargs) :: pt_string
        integer i

        nstrs=nargs
        pt_string=c_null_ptr
        do i = 1, nstrs
           array_strings(i) = trim(args(i))//c_null_char
           pt_string(i) = c_loc(array_strings(i))
        end do

        call ppidd_initialize_f2c(nstrs,pt_string)

      end subroutine ppidd_initialize
  
      end module ppidd_init_mod
