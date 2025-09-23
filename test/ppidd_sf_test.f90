! PPIDD SF standalone test suites

program main
use ppidd
implicit none
integer(c_int) helper_server_flag, num_per_helper, argc
type(c_ptr), target, dimension(1) :: argv

! Intitialize PPIDD library
argv(1) = c_null_ptr
argc = 0_c_int
call ppidd_initialize(argc,c_loc(argv),ppidd_impl_default,int(bit_size(0),c_int))
! one helper server on every node */
helper_server_flag=1_c_int
num_per_helper=0_c_int
call ppidd_helper_server(helper_server_flag,num_per_helper)
call ppidd_initialize_data()

call ppidd_sf_test()

! Terminate and Tidy up PPIDD
call ppidd_finalize()
end program main

subroutine ppidd_sf_test
use ppidd
implicit none
integer, parameter :: iout = 6

integer size
integer dimnsn
integer(c_int), parameter :: maxid=5 ! max number of outstanding I/O requests
parameter(dimnsn=8192)
double precision, target :: buffer(dimnsn,maxid)  ! need buffering for all maxid requests
double precision tt0,tt1,tcgtime
double precision, target :: ttw, ttr
!     integer stack, heap
integer(c_int) idlist(maxid)
character(len=80) errmsg
integer lerrmsg

integer nodeid, nnodes
integer i,start,end,j, chunk
integer(c_int) me, nproc, rc, curid, handle
character(8), parameter :: fname='ppidd_sf'

!...  Get the processes number and process id
nproc = ppidd_size()
me = ppidd_rank()
if (ppidd_impl_default.eq.ppidd_impl_no_mpi) then
  write(iout,*) 'Serial PPIDD, no tests performed.'
  return
end if
if(me.eq.0) write(iout,*) 'PPIDD initialized '
if(me.eq.0)then
   write(iout,*)
   write(iout,'(1x,2(a,i4))') 'Nprocs=',nproc,'   My proc=',me
   write(iout,*) 'Performing ppidd sf tests.'
   write(iout,*)
endif
!
curid = 0

size  = maxid*dimnsn*nproc
!
if(me.eq.0) write(iout,*) 'Creating shared file = ',fname
flush(6)
rc=ppidd_sf_create(fname//c_null_char,dble(16*size),dble(8*size),dble(8*dimnsn),handle)
!
call ppidd_barrier()
if(me.eq.0) write(iout,*) 'Writing and reading operations:'
flush(6)
call ppidd_barrier()
chunk = (size+nproc-1)/nproc
start = me*chunk+1
end = min((start+chunk-1),size)
tt0 = ppidd_wtime()

write(iout,*) 'me=',me,'writing:', start, end
!     everybody writes chunk of data
if(start.le.end) then
   do i = start, end,dimnsn
      do j = 1, min(dimnsn,(end-i+1))
         buffer(j,curid+1) = dble(i+j-1)
      enddo

      if(curid .eq. maxid)then
         rc=ppidd_sf_waitall(idlist,maxid)
         curid = 0
      endif
      curid = curid+1
      rc=ppidd_sf_write(handle,  8.0d0*dble(i-1), &
         8.0d0*dble(min(dimnsn,(end-i+1))), &
         c_loc(buffer(1,curid)), &
         idlist(curid))
      if (rc.ne.0)call ppidd_error('write failed'//c_null_char,rc)

   enddo
endif

rc=ppidd_sf_waitall(idlist,curid)
if(rc.ne.0)call ppidd_error('waitall failed'//c_null_char,rc)
curid = 0

tt1 = ppidd_wtime()
ttw = tt1 -tt0

call ppidd_gsum(1_c_int,c_loc(ttw),1_c_int,'max'//c_null_char)
call ppidd_barrier()
flush(6)


!     everybody reads different chunk of data
start = (nproc-me-1)*chunk+1
end = min((start+chunk-1),size)
write(iout,*) 'me=',me,'reading:', start, end
tt0 = ppidd_wtime()
do i = start,end,dimnsn

!           read and test data chunk by chunk
      rc=ppidd_sf_read(handle, 8.0d0*dble(i-1), &
         8.0d0*dble(min(dimnsn,(end-i+1))), c_loc(buffer), &
         idlist(1))
      if (rc.ne.0)then
         lerrmsg=len(errmsg)
         errmsg(lerrmsg:lerrmsg)=c_null_char
         call ppidd_sf_errmsg(rc,errmsg)
         write(iout,*) 'read at offset ',8.0d0*dble(i-1),' failed:',errmsg
         call ppidd_error('read failed'//c_null_char,rc)
      endif
      rc=ppidd_sf_wait(idlist(1))
      if (rc.ne.0)call ppidd_error('wait failed'//c_null_char,rc)

      do j = 1,min(dimnsn,(end-i+1))
         if(buffer(j,1).ne.dble(i+j-1)) then
            print *, me,buffer(j,1), i+j-1,i
            stop 'test failed'
         endif
      enddo
enddo
tt1 = ppidd_wtime()
ttr = tt1 -tt0

call ppidd_gsum(1_c_int,c_loc(ttr),1_c_int,'max'//c_null_char)
call ppidd_barrier()

flush(6)
if(me.eq.0)then
  write(iout,*)
  write(iout,*)'test passed ', 8*maxid*dimnsn,' bytes'
  write(iout,*) 8.0d-6*dble(maxid*dimnsn)/ttw,' MB/s write rate'
  write(iout,*) 8.0d-6*dble(maxid*dimnsn)/ttr,' MB/s read rate'
  write(iout,*)
  write(iout,*) 'Destroying shared file = ',fname
endif
rc=ppidd_sf_destroy(handle)

call ppidd_barrier()

if(me.eq.0) write(iout,*) 'PPIDD sf test successful.'
end subroutine ppidd_sf_test
