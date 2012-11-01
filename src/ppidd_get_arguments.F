#ifdef MOLPRO_f2003
#define PPIDD_f2003
#endif

      INTEGER FUNCTION ppidd_iargc()
      IMPLICIT NONE
#ifdef PPIDD_f2003
      INTEGER command_argument_count
      ppidd_iargc=command_argument_count()
#else
      INTEGER iargc
      ppidd_iargc=iargc()
#endif
      RETURN
      END FUNCTION ppidd_iargc

      SUBROUTINE ppidd_getarg(i,str)
      IMPLICIT NONE
      INTEGER          :: i
      CHARACTER(LEN=*) :: str
#if defined(__OPEN64__) || defined(__PATHSCALE__)
      INTEGER*4        :: j
#endif
#ifdef PPIDD_f2003
#if defined(__OPEN64__) || defined(__PATHSCALE__)
      j=i
      CALL get_command_argument(j,str)
#else
      CALL get_command_argument(i,str)
#endif
#else
      CALL getarg(i,str)
#endif
      RETURN
      END SUBROUTINE ppidd_getarg
