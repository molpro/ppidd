! PPIDD SF standalone test suites
c
c FNAME - filename for test program
c
#define FNAME   'ppidd_sf'

      program main
      implicit none
      logical helper_server_flag
      integer num_per_helper

! Intitialize PPIDD library
      call ppidd_initialize()
#ifdef NO_NXTVAL_SERVER
      helper_server_flag=.false.
      num_per_helper=-1
#else
! one helper server on every node */
      helper_server_flag=.true.
      num_per_helper=0
#endif
      call ppidd_helper_server(helper_server_flag,num_per_helper)
      call ppidd_initialize_data()

      call ppidd_sf_test()

! Terminate and Tidy up PPIDD
      call ppidd_finalize()
      end

      subroutine ppidd_sf_test
      implicit none
      integer iout
c
      integer size
      integer dimnsn
      integer maxid
      parameter(maxid=5)           ! max number of outstanding I/O requests
      parameter(dimnsn=8192)
      double precision buffer(dimnsn,maxid)  ! need buffering for all maxid requests
      double precision tt0,tt1,ttw, ttr, tcgtime
c     integer stack, heap
      integer idlist(maxid), curid
      character(len=80) errmsg
c
      integer nodeid, nnodes
      integer me,nproc,handle,i,start,end,j, rc, chunk

c...  Get the processes number and process id
      call ppidd_size(nproc)
      call ppidd_rank(me)
      iout=6
      if(me.eq.0) write(iout,*) 'PPIDD initialized '
      if(me.eq.0)then
         write(iout,*)
         write(iout,'(1x,2(a,i4))') 'Nprocs=',nproc,'   My proc=',me
         write(iout,*) 'Performing ppidd sf tests.'
         write(iout,*)
      endif
c
      curid = 0

      size  = maxid*dimnsn*nproc
c
      if(me.eq.0) write(iout,*) 'Creating shared file = ',
     >            trim(FNAME)
      call ppidd_flush6
      call ppidd_sf_create(FNAME,
     $     dble(16*size),dble(8*size),dble(8*dimnsn),handle,rc)
c
      call ppidd_barrier()
      if(me.eq.0) write(iout,*) 'Writing and reading operations:'
      call ppidd_flush6
      call ppidd_barrier()
      chunk = (size+nproc-1)/nproc
      start = me*chunk+1
      end = min((start+chunk-1),size)
      call ppidd_wtime(tt0)
c
      write(iout,*) 'me=',me,'writing:', start, end
c     everybody writes chunk of data
      if(start.le.end) then
         do i = start, end,dimnsn
            do j = 1, min(dimnsn,(end-i+1))
               buffer(j,curid+1) = dble(i+j-1)
            enddo

            if(curid .eq. maxid)then
               call ppidd_sf_waitall(idlist,maxid,rc)
               curid = 0
            endif
            curid = curid+1
            call ppidd_sf_write(handle,  8.0d0*dble(i-1),
     &           8.0d0*dble(min(dimnsn,(end-i+1))), buffer(1,curid),
     &           idlist(curid),rc)
            if (rc.ne.0)call ppidd_error('write failed',rc)

         enddo
      endif
c
      call ppidd_sf_waitall(idlist,curid,rc)
      if(rc.ne.0)call ppidd_error('waitall failed',rc)
      curid = 0

      call ppidd_wtime(tt1)
      ttw = tt1 -tt0
c
      call ppidd_gsum(1,ttw,1,'max')
      call ppidd_barrier()
      call ppidd_flush6
c
c
c     everybody reads different chunk of data
      start = (nproc-me-1)*chunk+1
      end = min((start+chunk-1),size)
      write(iout,*) 'me=',me,'reading:', start, end
      call ppidd_wtime(tt0)
      do i = start,end,dimnsn

c           read and test data chunk by chunk
            call ppidd_sf_read(handle, 8.0d0*dble(i-1),
     .          8.0d0*dble(min(dimnsn,(end-i+1))), buffer, idlist(1),rc)
            if (rc.ne.0)then
               call ppidd_sf_errmsg(rc,errmsg)
               write(iout,*) 'read at offset ',8.0d0*dble(i-1),
     .                       ' failed:',errmsg
               call ppidd_error('read failed',rc)
            endif
            call ppidd_sf_wait(idlist(1),rc)
            if (rc.ne.0)call ppidd_error('wait failed',rc)
c
            do j = 1,min(dimnsn,(end-i+1))
               if(buffer(j,1).ne.dble(i+j-1)) then
                  print *, me,buffer(j,1), i+j-1,i
                  stop 'test failed'
               endif
            enddo
      enddo
      call ppidd_wtime(tt1)
      ttr = tt1 -tt0
c
      call ppidd_gsum(1,ttr,1,'max')
      call ppidd_barrier()
c
      call ppidd_flush6
      if(me.eq.0)then
        write(iout,*)
        write(iout,*)'test passed ', 8*maxid*dimnsn,' bytes'
        write(iout,*) 8.0d-6*dble(maxid*dimnsn)/ttw,' MB/s write rate'
        write(iout,*) 8.0d-6*dble(maxid*dimnsn)/ttr,' MB/s read rate'
        write(iout,*)
        write(iout,*) 'Destroying shared file = ',trim(FNAME)
      endif
      call ppidd_sf_destroy(handle,rc)
c
      call ppidd_barrier()

      if(me.eq.0) write(iout,*) 'PPIDD sf test successful.'
c
      return
      end
