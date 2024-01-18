! PPIDD standalone test suites

      program main
      use ppidd
      implicit none
      integer(c_int)     :: helper_server_flag
      integer(c_int64_t) :: num_per_helper
      integer(c_int64_t) :: iprocs
      integer(c_int)     :: argc=0_c_int
      type(c_ptr), target, dimension(1) :: argv
      integer iout
#ifdef PPIDD_USE_MA
      integer ma_heap, idtype
      logical status
      integer, parameter :: heap=400*400*4
      integer, parameter :: fudge=100
      integer, parameter :: stack=400*400*4
#endif

!  Intitialize PPIDD library
!  Intitialize a message passing library
!  Initialize GA
      argv(1)=c_null_ptr
      call ppidd_initialize(argc,c_loc(argv),ppidd_impl_default)
#ifdef NO_NXTVAL_SERVER
      helper_server_flag=0_c_int
      num_per_helper=-1_c_int
#else
! one helper server on every node */
      helper_server_flag=1_c_int64_t
      num_per_helper=0_c_int64_t
#endif
      call ppidd_helper_server(helper_server_flag,num_per_helper)
      call ppidd_initialize_data()
c
      call ppidd_rank(iprocs)
      iout=6
      if(iprocs.eq.0)then
         write(iout,*)
#ifdef MPI2
         write(iout,*) 'PPIDD initialized (MPI-2 version of PPIDD).'
#else
         write(iout,*) 'PPIDD initialized (GA version of PPIDD).'
#endif
         call ppidd_flush6
      endif
      call ppidd_barrier()

#if defined(PPIDD_USE_MA)
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      ma_heap=heap
      call ppidd_uses_ma(status)
      IF (status) then
        ma_heap=heap+fudge
        if(iprocs.eq.0) write(iout,11) ma_heap
11      format(' ppidd_uses_ma= true, calling ppidd_ma_init with',
     >         ' heap size=',i15)
      else
        if(iprocs.eq.0) write(iout,22)
22      format(' ppidd_uses_ma=false, calling ppidd_ma_init with',
     >         ' nominal heap.')
      END IF
c... idtype=1: fortran double precision type
      idtype=1
      call ppidd_ma_init(idtype,stack,ma_heap,status)
      if (.not. status)
     > call ppidd_error('ppidd_ma_init failed'//c_null_char,
     >                  -1_c_int)

      if(iprocs.eq.0) write(iout,*) 'MA initialized.'
#endif

! Call ppiddtest
      call ppiddtest(helper_server_flag.eq.1)

! Terminate and Tidy up PPIDD
      call ppidd_finalize()
      end

      subroutine ppiddtest(helper_server_flag)
      use ppidd
      implicit double precision (a-h,o-z)
      integer(c_int64_t) :: iprocs, nproc, nproc_total
      logical helper_server_flag

c...  Get the processes number and process id
      call ppidd_size(nproc)
      call ppidd_rank(iprocs)
      call ppidd_size_all(nproc_total)
c
      iout=6
c
c... Print messages
c
      if(iprocs.eq.0)then
         write(iout,*)
         write(iout,19) nproc_total,nproc,iprocs
         write(iout,*) 'Performing ppidd tests.'
19       format(1x,'Total procs=',i4,',  Nprocs(compute processes)=',
     >      i4,',  My proc=',i4)
         write(iout,29)
29       format(1x,'SUCCESSFUL TEST SHOULD END WITH ',
     >      ''' All PPIDD tests successful.''')
         write(iout,*)
         call ppidd_flush6
      endif

! Synchronize processes (a barrier) before specific tests
#if defined(MPI2) || defined(GA_MPI)
      call ppidd_barrier()
c... Checking if C data types match Fortran
      if(iprocs.eq.0)then
        call ppidd_checkdtype
      endif
      call ppidd_barrier()
c... Performing PPIDD basic operation test
      call get_current_times(cpua,wtimea)
      call ppidd_basictest(helper_server_flag)
      call ppidd_barrier()
      call get_current_times(cpub,wtimeb)
      if(iprocs.eq.0)write(iout,*)
      if(iprocs.eq.0)write(iout,100) 'ppidd_basictest:',
     >   cpub-cpua,wtimeb-wtimea
100   format(1x,'Time spent on ',a,t49,' cpu=',f16.4,' sec,',
     > '  wall time=',f16.4,' sec')

      call ppidd_barrier()
c... Performing PPIDD mutex test
      call get_current_times(cpua,wtimea)
      call ppidd_mutex_test(helper_server_flag)
      call ppidd_barrier()
      call get_current_times(cpub,wtimeb)
      if(iprocs.eq.0)write(iout,*)
      if(iprocs.eq.0)write(iout,100) 'ppidd_mutex_test:',
     >   cpub-cpua,wtimeb-wtimea
      call ppidd_barrier()
c... Performing PPIDD shared counter test
      call get_current_times(cpua,wtimea)
      if (helper_server_flag) then
        call ppidd_sharedcounter_test
      else
        if(iprocs.eq.0) write(iout,*)
     >    'Helper server is disabled. Skip the shared counter test.'
      end if
      call ppidd_barrier()
      call get_current_times(cpub,wtimeb)
      if(iprocs.eq.0)write(iout,*)
      if(iprocs.eq.0)write(iout,100) 'ppidd_sharedcounter_test:',
     >   cpub-cpua,wtimeb-wtimea
c
      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,*) 'All PPIDD tests successful.'
#else
      write(iout,*) 'All PPIDD tests are skipped in serial case.'
#endif
c
      return
      end


! Checking whether C types match Fortran types
      subroutine ppidd_checkdtype
      use ppidd
      implicit double precision (a-h,o-z)
      logical ok_int,ok_log,ok_dbl

      interface
       subroutine sizeofctype(sizeint,sizelog,sizedouble,sizefloat)
     >            bind(C,name='sizeofctype')
       use iso_c_binding, only : c_double
       real(c_double), intent(out) :: sizeint
       real(c_double), intent(out) :: sizelog
       real(c_double), intent(out) :: sizedouble
       real(c_double), intent(out) :: sizefloat
       end subroutine sizeofctype
      end interface

      write(6,*) 'Checking whether C types match Fortran types:'
      call sizeofctype(dbyte_int1,dbyte_log1,dbyte_dbl1,dbyte_float)
      nbyte_int1=nint(dbyte_int1)
      nbyte_log1=nint(dbyte_log1)
      nbyte_dbl1=nint(dbyte_dbl1)
      nbyte_float=nint(dbyte_float)
      write(6,*) 'The size (in bytes) of some C data types:'
      write(6,10) 'float:', nbyte_float
      write(6,10) 'double:', nbyte_dbl1
      write(6,10) 'int64_t:', nbyte_int1
      write(6,10) 'fortlogical:', nbyte_log1
10    format(1x,'Size of C ',a,t35,i4,' bytes')
      call sizeoffortype(nbyte_int2,nbyte_log2,nbyte_dbl2,nbyte_real)
      write(6,*) 'The size (in bytes) of some Fortran data types:'
      write(6,20) 'Real:', nbyte_real
      write(6,20) 'Doule Precision:', nbyte_dbl2
      write(6,20) 'Integer:', nbyte_int2
      write(6,20) 'Logical:', nbyte_log2
20    format(1x,'Size of Fortran ',a,t35,i4,' bytes')
      ok_int=nbyte_int1.eq.nbyte_int2
      ok_log=nbyte_log1.eq.nbyte_log2
      ok_dbl=nbyte_dbl1.eq.nbyte_dbl2
      if( ok_int .and. ok_log .and. ok_dbl) then
        write(6,30)
30      format(1x,'All data types (Integer, Logical, and ',
     >         'Double Precision) match between C and Fortran.')
      else
        if (.not.ok_int) write(6,31) 'Integer'
        if (.not.ok_log) write(6,31) 'Logical'
        if (.not.ok_dbl) write(6,31) 'Double Precision'
31      format(1x,'ERROR: Data type ( ',a,t38,') does not match between'
     >           ,' C and Fortran types.')
        write(6,32)
32      format(1x,'Please check INTEGER (related to CFLAGS/FCFLAGS)',
     >            ' option, and rebuild the library.')
        call ppidd_error('Data types failed to match'//c_null_char,
     >                   0_c_int)
      endif

      end


      subroutine ppidd_basictest(helper_server_flag)
      use ppidd
      implicit double precision (a-h,o-z)
      logical helper_server_flag
      double precision, allocatable, target :: buff(:)
      double precision totsize,bufsize
      double precision, target :: z1
      double precision cpu1,cpu2,cpu,speed
      integer lentot,lenbuf,len,ioff,ilo,ihi
      integer lenseg,lenlat
      integer(c_int64_t) :: nprocs, iprocs
      integer maxprc,nloop,iout,i,ipr
      integer proc_ltop,proc_rbot,dest_src,sync,lenre,srcre
      integer(c_int) dtype
      integer storetype
      integer ihandle,ihandlat
      character(len=80) :: name
      integer(c_int) :: ok
      logical verbose
      logical arcca_test_flag

      verbose=.false.
      totsize=0.0d0
      bufsize=0.0d0
      maxprc=0
      nloop=10

      iout=6
! ARCCA Test Flag
      arcca_test_flag=.false.
!      arcca_test_flag=.true.

! ARCCA test?
      if (arcca_test_flag) then
        nrep_1side=200000
        nrep_send =2000000
        nrep_bcast=1000000
        nrep_gsum =2000000
      else
        nrep_1side=20000
        nrep_send =200000
        nrep_bcast=100000
        nrep_gsum =200000
      end if

      call ppidd_size(nprocs)
      call ppidd_rank(iprocs)
      if(maxprc.eq.0) maxprc=nprocs
      maxprc=min(maxprc,nprocs)
      if(nloop.eq.0) nloop=1
      if(dabs(totsize).lt.1.0d0) totsize=32d0
      if(arcca_test_flag) totsize=128d0  !ARCCA
      if(dabs(bufsize).lt.1.0d0) bufsize=16d0
      lentot=nint(totsize)*1024*128
      lenbuf=nint(bufsize)*1024*128
      if(iprocs.eq.0) write(iout,1) totsize,lentot,bufsize,lenbuf
1     format(/' Starting PPIDD basic test with GA-size=',f8.2,
     >        ' MB (',i9,' words),  Buffer size=',f8.2,
     >        ' MB (',i9,' words)'/)
c
c... initialise the buffer space
      ALLOCATE( buff(lenbuf) )
      do i=1,lenbuf
        buff(i)=10.01d0+dble(i)*0.001d0
      end do

#ifdef MPI2
      if (helper_server_flag) then
        nloop_array=2
      else
        nloop_array=1
      end if
#else
      nloop_array=1
#endif
c     dtype=1_c_int: double precision
      do loop_array=1,nloop_array
      if (nloop_array.eq.1) then
        if(iprocs.eq.0) write(iout,101) ' '
101     format(/,a,'Array data are stored across distributed compute'
     >  ,' processes:')
        dtype=1_c_int
        storetype=0
      else
c... case (1)
        if (loop_array.eq.1) then
          if(iprocs.eq.0) write(iout,101) ' (1) '
          dtype=1_c_int
          storetype=0
c... case (2)
        else
          if(iprocs.eq.0) write(iout,102)
102       format(/' (2) Array data are stored on helper processes:')
          dtype=1_c_int
          storetype=1
        end if
      end if
!      if (loop_array.eq.1 .and. nloop_array.eq.2) cycle  !skip case (1) for MPI2
!      if (loop_array.eq.2 .and. nloop_array.eq.2) cycle  !skip case (2) for MPI2
c
c... ppidd_create
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      ok=ppidd_create('ppiddarray'//c_null_char,lentot,dtype,
     >                storetype,ihandle)
      call get_current_times(cpu2,wtime2)
      cpu=cpu2-cpu1
      wtime=wtime2-wtime1
      call ppidd_inquire_name(ihandle,name)
      if(iprocs.eq.0)write(iout,*)
     >  'ppidd has created global array: ',name(1:len_trim(name))
      if(iprocs.eq.0)write(iout,10) 'ppidd_create:',cpu,wtime
10    format(1x,a,t20,'cpu=',f9.4,' sec,',
     >       '  wall time=',f9.4,' sec',:,
     >       ',  speed=', f12.2,' MB/sec')
      call ppidd_flush6
cc    create global data struct for testing latency, lenseg elements per proc evenly
      lenseg=1
      lenlat=nprocs*lenseg
      ok=ppidd_create('ppiddlat'//c_null_char,lenlat,dtype,storetype,
     >                ihandlat)
c
c... ppidd_zero
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      do i=1,nloop
       ok=ppidd_zero(ihandle)
      end do
      call get_current_times(cpu2,wtime2)
      cpu=cpu2-cpu1
      wtime=wtime2-wtime1
      if(iprocs.eq.0) write(iout,10) 'ppidd_zero:',cpu,wtime
      call ppidd_flush6
      ok=ppidd_zero(ihandlat)
c
c... ppidd_put
c... measure the latency
c  (1)  For integer global data structure stored on helper process, measure the time for operating 2 elements
c    since 1 element operation is a special case
c  (2)  But for double precision global data structure stored on helper process, measure the time for operating 1 element
      nrep=nrep_1side/(nprocs*nprocs)
      if(nrep.eq.0) nrep=1
c ppidd_put latency trial test, in order to eliminate unstable behavior
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            do inum=1,nprocs
              ok=ppidd_put(ihandlat,inum,inum,c_loc(buff(inum)))
            end do
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=wtimes/dble(nrep*nprocs*nprocs)
      if(iprocs.eq.0) write(iout,49) 'PPIDD_PUT=',wtimes,wlat*1.0d6
49    format(1x,'TIME FOR LATENCY TEST ',a,t41,f10.2,' SEC,  LATENCY=',
     >          f10.2,' MICROSEC (trial test)')

c ppidd_put latency standard test
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            do inum=1,nprocs
              ok=ppidd_put(ihandlat,inum,inum,c_loc(buff(inum)))
            end do
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=wtimes/dble(nrep*nprocs*nprocs)
      if(iprocs.eq.0) write(iout,50) 'PPIDD_PUT=',wtimes,wlat*1.0d6
50    format(1x,'TIME FOR LATENCY TEST ',a,t41,f10.2,' SEC,  LATENCY=',
     >          f10.2,' MICROSEC')


c ppidd_put bandwidth trial test, in order to eliminate unstable behavior
      nloop_save=nloop
      nloop=1
      cpus=0.0
      wtimes=0.0
      nopr=0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        call get_current_times(cpu1,wtime1)
        if(iprocs.eq.ipr) then
          do i=1,nloop
          nrest=lentot
          do ioff=0,lentot-1,lenbuf
            len=min(lenbuf,nrest)
            ilo=ioff+1
            ihi=ioff+len
            ok=ppidd_put(ihandle,ilo,ihi,c_loc(buff))
            nopr=nopr+1
            nrest=nrest-len
          end do
          end do
        end if
        call ppidd_barrier()
        call get_current_times(cpu2,wtime2)
        cpu=cpu2-cpu1
        wtime=wtime2-wtime1
        cpus=cpus+cpu
        wtimes=wtimes+wtime
      end do
      speed=0.0d0
      if(wtimes.gt.0d0) speed=dble(maxprc*nloop)*totsize/wtimes
      if(iprocs.eq.0)
     >   write(iout,10) 'ppidd_put(trial):',cpus,wtimes,speed
      call ppidd_flush6
      call ppidd_barrier()
      nloop=nloop_save


c... measure the bandwidth
      cpus=0.0
      wtimes=0.0
      nopr=0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        call get_current_times(cpu1,wtime1)
        if(iprocs.eq.ipr) then
          do i=1,nloop
          nrest=lentot
          do ioff=0,lentot-1,lenbuf
            len=min(lenbuf,nrest)
            ilo=ioff+1
            ihi=ioff+len
            ok=ppidd_put(ihandle,ilo,ihi,c_loc(buff))
            nopr=nopr+1
            nrest=nrest-len
          end do
          end do
        end if
        call ppidd_barrier()
        call get_current_times(cpu2,wtime2)
        cpu=cpu2-cpu1
        wtime=wtime2-wtime1
        cpus=cpus+cpu
        wtimes=wtimes+wtime
        if(iprocs.eq.0 .and. verbose) then
          speed=0.0d0
          if(wtime.gt.0d0) speed=totsize*dble(nloop)/wtime
          write(iout,22) 'ppidd_put:',ipr,cpu,wtime,speed
22        format(1x,a,t15,'  iprocs=',i3,'  cpu=',f9.4,' sec,',
     >         '  wall time=',f9.4,' sec,  speed=',
     >         f8.2,' MB/sec')
          call ppidd_flush6
        end if
      end do
      speed=0.0d0
      if(wtimes.gt.0d0) then
        speed=dble(maxprc*nloop)*totsize/wtimes
        bdw=dble(maxprc*nloop)*totsize/(wtimes-wlat*dble(nopr*nprocs))
      end if
      if(iprocs.eq.0.and.verbose) write(iout,23) 'ppidd_put:',speed
23    format(1x,a,t18,'  average speed=',f8.2,' MB/sec,')
      if(iprocs.eq.0) then
        write(iout,55) 'ppidd_put:',wlat*dble(nopr*nprocs),nopr
        write(iout,10) 'ppidd_put:',cpus,wtimes,speed
        write(iout,60) 'ppidd_put:',cpus,wtimes,wlat*1.0d6,bdw
      end if
55    format(1x,a,t20,' total latency time=',f9.4,' sec,  nopr=',i6)
60    format(1x,a,t20,'cpu=',f9.4,' sec,',
     >       '  wall time=',f9.4,' sec',:,
     >     ',  latency=',f10.2,' microsec, bandwidth=', f12.2,' MB/sec')
      call ppidd_flush6

c
c... ppidd_get
c... measure the latency
c ppidd_get latency trial test, in order to eliminate unstable behavior
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            do inum=1,nprocs
              ok=ppidd_get(ihandlat,inum,inum,c_loc(buff(inum)))
            end do
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=wtimes/dble(nrep*nprocs*nprocs)
      if(iprocs.eq.0) write(iout,49) 'PPIDD_GET=',wtimes,wlat*1.0d6

c ppidd_get latency standard test
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            do inum=1,nprocs
              ok=ppidd_get(ihandlat,inum,inum,c_loc(buff(inum)))
            end do
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=wtimes/dble(nrep*nprocs*nprocs)
      if(iprocs.eq.0) write(iout,50) 'PPIDD_GET=',wtimes,wlat*1.0d6

c... measure the bandwidth
      cpus=0.0
      wtimes=0.0
      nopr=0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        call get_current_times(cpu1,wtime1)
        if(iprocs.eq.ipr) then
          do i=1,nloop
          nrest=lentot
          do ioff=0,lentot-1,lenbuf
            len=min(lenbuf,nrest)
            ilo=ioff+1
            ihi=ioff+len
            ok=ppidd_get(ihandle,ilo,ihi,c_loc(buff))
            nopr=nopr+1
            nrest=nrest-len
          end do
          end do
        end if
        call ppidd_barrier()
        call get_current_times(cpu2,wtime2)
        cpu=cpu2-cpu1
        wtime=wtime2-wtime1
        cpus=cpus+cpu
        wtimes=wtimes+wtime
        if(iprocs.eq.0 .and. verbose) then
          speed=0.0d0
          if(wtime.gt.0d0) speed=totsize*dble(nloop)/wtime
          write(iout,22) 'ppidd_get:',ipr,cpu,wtime,speed
          call ppidd_flush6
        end if
      end do
      speed=0.0d0
      if(wtimes.gt.0d0) then
        speed=dble(maxprc*nloop)*totsize/wtimes
        bdw=dble(maxprc*nloop)*totsize/(wtimes-wlat*dble(nopr*nprocs))
      end if
      if(iprocs.eq.0.and.verbose) write(iout,23) 'ppidd_get:',speed
      if(iprocs.eq.0) then
        write(iout,55) 'ppidd_get:',wlat*dble(nopr*nprocs),nopr
        write(iout,10) 'ppidd_get:',cpus,wtimes,speed
        write(iout,60) 'ppidd_get:',cpus,wtimes,wlat*1.0d6,bdw
      end if
      call ppidd_flush6
c
c... ppidd_acc
c... measure the latency
      iz1=1
      z1=1.0d0
c ppidd_acc latency trial test, in order to eliminate unstable behavior
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            do inum=1,nprocs
              ok=ppidd_acc(ihandlat,inum,inum,c_loc(buff(inum)),
     >                     c_loc(z1))
            end do
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=wtimes/dble(nrep*nprocs*nprocs)
      if(iprocs.eq.0) write(iout,49) 'PPIDD_ACC=',wtimes,wlat*1.0d6

c ppidd_acc latency standard test
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            do inum=1,nprocs
              ok=ppidd_acc(ihandlat,inum,inum,c_loc(buff(inum)),
     >                     c_loc(z1))
            end do
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=wtimes/dble(nrep*nprocs*nprocs)
      if(iprocs.eq.0) write(iout,50) 'PPIDD_ACC=',wtimes,wlat*1.0d6

c... measure the bandwidth
      cpus=0.0
      wtimes=0.0
      nopr=0
      do ipr=0,maxprc-1
        call ppidd_barrier()
        call get_current_times(cpu1,wtime1)
        if(iprocs.eq.ipr) then
          do i=1,nloop
          nrest=lentot
          do ioff=0,lentot-1,lenbuf
            len=min(lenbuf,nrest)
            ilo=ioff+1
            ihi=ioff+len
            ok=ppidd_acc(ihandle,ilo,ihi,c_loc(buff),c_loc(z1))
            nopr=nopr+1
            nrest=nrest-len
          end do
          end do
        end if
        call ppidd_barrier()
        call get_current_times(cpu2,wtime2)
        cpu=cpu2-cpu1
        wtime=wtime2-wtime1
        cpus=cpus+cpu
        wtimes=wtimes+wtime
        if(iprocs.eq.0 .and. verbose) then
          speed=0.0d0
          if(wtime.gt.0d0) speed=totsize*dble(nloop)/wtime
          write(iout,22) 'ppidd_acc:',ipr,cpu,wtime,speed
          call ppidd_flush6
        end if
      end do
      speed=0.0
      if(wtimes.gt.0d0) then
        speed=dble(maxprc*nloop)*totsize/wtimes
        bdw=dble(maxprc*nloop)*totsize/(wtimes-wlat*dble(nopr*nprocs))
      end if
      if(iprocs.eq.0.and.verbose) write(iout,23) 'ppidd_acc:',speed
      if(iprocs.eq.0) then
        write(iout,55) 'ppidd_acc:',wlat*dble(nopr*nprocs),nopr
        write(iout,10) 'ppidd_acc:',cpus,wtimes,speed
        write(iout,60) 'ppidd_acc:',cpus,wtimes,wlat*1.0d6,bdw
      end if
      call ppidd_flush6
c
c... ppidd_destroy
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      ok=ppidd_destroy(ihandle)
      call get_current_times(cpu2,wtime2)
      cpu=cpu2-cpu1
      wtime=wtime2-wtime1
      if(iprocs.eq.0) write(iout,10) 'ppidd_destroy:',cpu,wtime
      call ppidd_flush6

      ok=ppidd_destroy(ihandlat)

      end do
c...  end of loop_array

c
c... ppidd_send and ppidd_recv
      if(iprocs.eq.0)write(iout,*)
      proc_ltop=maxprc/2
      proc_rbot=(maxprc-1)/2
      dest_src=maxprc-1-iprocs
c     data type: double precision; sync: synchronous
      dtype=1_c_int
      sync=1

c.... measure latency
      n=1
      nrep=nrep_send*2/nprocs
      if(nrep.eq.0) nrep=1
c ppidd_send/recv latency trial test, in order to eliminate unstable behavior
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,proc_ltop-1
        ipr_pair=maxprc-1-ipr
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            call ppidd_send(c_loc(buff(iprocs)),n,dtype,dest_src,sync)
          end do
        end if
        if(iprocs.eq.ipr_pair) then
          do irep=1,nrep
            call ppidd_recv(c_loc(buff(iprocs)),n,dtype,dest_src,lenre,
     >                      srcre,sync)
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=0.0
      if (proc_ltop.ne.0) wlat=wtimes/dble(nrep*proc_ltop)
      if(iprocs.eq.0)write(iout,49) 'PPIDD_SEND/RECV=',wtimes,wlat*1.d6

c ppidd_send/recv latency standard test
      call ppidd_barrier()
      wtimes=0.0
      do ipr=0,proc_ltop-1
        ipr_pair=maxprc-1-ipr
        call ppidd_barrier()
        wtime1=ppidd_wtime()
        if(iprocs.eq.ipr) then
          do irep=1,nrep
            call ppidd_send(c_loc(buff(iprocs)),n,dtype,dest_src,sync)
          end do
        end if
        if(iprocs.eq.ipr_pair) then
          do irep=1,nrep
            call ppidd_recv(c_loc(buff(iprocs)),n,dtype,dest_src,lenre,
     >                      srcre,sync)
          end do
        end if
        call ppidd_barrier()
        wtime2=ppidd_wtime()
        wtimes=wtimes+wtime2-wtime1
      end do
      call ppidd_barrier()
      wlat=0.0
      if (proc_ltop.ne.0) wlat=wtimes/dble(nrep*proc_ltop)
      if(iprocs.eq.0)write(iout,50) 'PPIDD_SEND/RECV=',wtimes,wlat*1.d6

c.... measure bandwidth
      cpus=0.0
      wtimes=0.0
      dtype=1_c_int
      nopr=0
!      nloop_save=nloop
!      nloop=nloop*8
      do ipr=0,proc_ltop-1
        ipr_pair=maxprc-1-ipr
        call ppidd_barrier()
        call get_current_times(cpu1,wtime1)
        if(iprocs.eq.ipr) then
          do i=1,nloop
          nrest=lentot
          do ioff=0,lentot-1,lenbuf
            len=min(lenbuf,nrest)
            call ppidd_send(c_loc(buff),len,dtype,dest_src,sync)
            nopr=nopr+1
            nrest=nrest-len
          end do
          end do
        end if
        if(iprocs.eq.ipr_pair) then
          do i=1,nloop
          nrest=lentot
          do ioff=0,lentot-1,lenbuf
            len=min(lenbuf,nrest)
            call ppidd_recv(c_loc(buff),len,dtype,dest_src,lenre,srcre,
     >                      sync)
            nrest=nrest-len
          end do
          end do
        end if
        call ppidd_barrier()
        call get_current_times(cpu2,wtime2)
        cpu=cpu2-cpu1
        wtime=wtime2-wtime1
        cpus=cpus+cpu
        wtimes=wtimes+wtime
        if(iprocs.eq.0 .and. verbose) then
          speed=0.0d0
          if(wtime.gt.0d0) speed=totsize*dble(nloop)/wtime
          write(iout,22) 'ppidd_send/recv:',ipr,cpu,wtime,speed
          call ppidd_flush6
        end if
      end do
      speed=0.0
      if(wtimes.gt.0d0) then
         speed=totsize*dble(proc_ltop*nloop)/wtimes
         bdw=totsize*dble(proc_ltop*nloop)/(wtimes-
     >       wlat*dble(nopr*proc_ltop))
      end if
      if(iprocs.eq.0.and.verbose)
     >  write(iout,23) 'ppidd_send/recv:',speed
      if(iprocs.eq.0) write(iout,55) 'ppidd_send/recv:',
     >  wlat*dble(nopr*proc_ltop),nopr
      if(iprocs.eq.0) write(iout,10) 'ppidd_send/recv:',
     >  cpus,wtimes,speed
      if(iprocs.eq.0) write(iout,60) 'ppidd_send/recv:',
     >  cpus,wtimes,wlat*1.0d6,bdw
      call ppidd_flush6
!      nloop=nloop_save
c
c... ppidd_brdcst
      n=1
      nrep=nrep_bcast/nprocs
      if(nrep.eq.0) nrep=1
c ppidd_bcast latency trial test, in order to eliminate unstable behavior
      call ppidd_barrier()
      wtimes=0.0
      call ppidd_barrier()
      wtime1=ppidd_wtime()
      do irep=1,nrep
        do ipr=0,nprocs-1
          call ppidd_bcast(c_loc(buff),n,1_c_int,ipr)
        end do
      end do
      call ppidd_barrier()
      wtime2=ppidd_wtime()
      wtimes=wtime2-wtime1
      wlat=wtimes/dble(nrep*nprocs)
      if(iprocs.eq.0) write(iout,49) 'PPIDD_BRDCST=',wtimes,wlat*1.0d6

c ppidd_bcast latency standard test
      call ppidd_barrier()
      wtimes=0.0
      call ppidd_barrier()
      wtime1=ppidd_wtime()
      do irep=1,nrep
        do ipr=0,nprocs-1
          call ppidd_bcast(c_loc(buff),n,1_c_int,ipr)
        end do
      end do
      call ppidd_barrier()
      wtime2=ppidd_wtime()
      wtimes=wtime2-wtime1
      wlat=wtimes/dble(nrep*nprocs)
      if(iprocs.eq.0) write(iout,50) 'PPIDD_BRDCST=',wtimes,wlat*1.0d6

c ppidd_bcast bandwidth test
      cpus=0.0
      wtimes=0.0
      nopr=0
      totsize_save=totsize
      lentot_save=lentot
      totsize=totsize*8.0d0/dble(nprocs)
      lentot=nint(totsize)*1024*128
      do ipr=0,maxprc-1
        call ppidd_barrier()
        call get_current_times(cpu1,wtime1)
        do i=1,nloop
        nrest=lentot
        do ioff=0,lentot-1,lenbuf
          len=min(lenbuf,nrest)
          call ppidd_bcast(c_loc(buff),len,1_c_int,ipr)
          nopr=nopr+1
          nrest=nrest-len
        end do
        end do
        call ppidd_barrier()
        call get_current_times(cpu2,wtime2)
        cpu=cpu2-cpu1
        wtime=wtime2-wtime1
        cpus=cpus+cpu
        wtimes=wtimes+wtime
        if(iprocs.eq.0 .and. verbose) then
          speed=0.0d0
          if(wtime.gt.0d0) speed=totsize*dble(nloop)/wtime
          write(iout,22) 'ppidd_brdcst:',ipr,cpu,wtime,speed
          call ppidd_flush6
        end if
      end do
      speed=0.0
      bdw=0.0d0
      if(wtimes.gt.0d0) speed=dble(maxprc*nloop)*totsize/wtimes
      if(wtimes-wlat*dble(nopr) .gt. 0d0) then
        bdw=dble(maxprc*nloop)*totsize/(wtimes-wlat*dble(nopr))
      end if
      if(iprocs.eq.0.and.verbose)
     >    write(iout,23) 'ppidd_brdcst:',speed
      if(iprocs.eq.0) then
        write(iout,55) 'ppidd_brdcst:',wlat*dble(nopr),nopr
        write(iout,10) 'ppidd_brdcst:',cpus,wtimes,speed
        write(iout,60) 'ppidd_brdcst:',cpus,wtimes,wlat*1.0d6,bdw
      end if
      call ppidd_flush6
      totsize=totsize_save
      lentot=lentot_save
c
c... ppidd_gsum
      n=1
      nrep=nrep_gsum/nprocs
      if(nrep.eq.0) nrep=1
c ppidd_gsum latency trial test, in order to eliminate unstable behavior
      call ppidd_barrier()
      wtimes=0.0
      call ppidd_barrier()
      wtime1=ppidd_wtime()
      do irep=1,nrep
        call ppidd_gsum(1_c_int,c_loc(buff),n,'+'//c_null_char)
      end do
      call ppidd_barrier()
      wtime2=ppidd_wtime()
      wtimes=wtime2-wtime1
      wlat=wtimes/dble(nrep)
      if(iprocs.eq.0) write(iout,49) 'PPIDD_GSUM=',wtimes,wlat*1.0d6

c ppidd_gsum latency standard test
      call ppidd_barrier()
      wtimes=0.0
      call ppidd_barrier()
      wtime1=ppidd_wtime()
      do irep=1,nrep
        call ppidd_gsum(1_c_int,c_loc(buff),n,'+'//c_null_char)
      end do
      call ppidd_barrier()
      wtime2=ppidd_wtime()
      wtimes=wtime2-wtime1
      wlat=wtimes/dble(nrep)
      if(iprocs.eq.0) write(iout,50) 'PPIDD_GSUM=',wtimes,wlat*1.0d6

c ppidd_gsum bandwidth standard test
      nopr=0
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      do i=1,nloop
      nrest=lentot
      do ioff=0,lentot-1,lenbuf
        len=min(lenbuf,nrest)
        call ppidd_gsum(1_c_int,c_loc(buff),len,'+'//c_null_char)
        nopr=nopr+1
        nrest=nrest-len
      end do
      end do
      call ppidd_barrier()
      call get_current_times(cpu2,wtime2)
      cpus=cpu2-cpu1
      wtimes=wtime2-wtime1
      speed=0.0
      if(wtimes.gt.0d0) then
        speed=totsize*dble(nloop)/wtimes
        bdw=totsize*dble(nloop)/(wtimes-wlat*dble(nopr))
      end if
      if(iprocs.eq.0.and.verbose) write(iout,23) 'ppidd_gsum:',speed
      if(iprocs.eq.0) write(iout,55) 'ppidd_gsum:',wlat*dble(nopr),nopr
      if(iprocs.eq.0)write(iout,10) 'ppidd_gsum:',cpus,wtimes,speed
      if(iprocs.eq.0) write(iout,60) 'ppidd_gsum:',
     >  cpus,wtimes,wlat*1.0d6,bdw
      call ppidd_flush6
c
c... release the memory in buffer space
      IF( ALLOCATED(buff) ) DEALLOCATE( buff )

      return
      end

c... get current CPU time from Fortran intrinsic function and wall time from MPI2/GA function
      subroutine get_current_times(cpu,wtime)
      use ppidd
      implicit double precision (a-h,o-z)
      double precision cpu,wtime

      cpu=dble(second())
      wtime=ppidd_wtime()

      return
      end



      subroutine ppidd_sharedcounter_test
      use ppidd
      implicit double precision (a-h,o-z)
      integer,          allocatable, target :: task_array(:)
      double precision, allocatable, target :: sum_array(:)
      double precision, allocatable, target :: tcpu_array1(:)
      double precision, allocatable, target :: twt_array1(:)
      double precision, allocatable, target :: twt_array2(:)
      double precision cpu1,cpu2
      double precision wtime1,wtime2
      double precision sumtemp,diffm
      integer(c_int64_t) :: nprocs,iprocs
      integer iout,i
      integer num_tasks,nscale_task,npart
      integer nval,junk
      logical verbose
      logical arcca_test_flag
c
c... set task number and task scale
c    please change these two numbers if carrying out large-scale, many-task test
c    eg. num_tasks=1000 nscale_task=1000 (Please be aware that it may take hours depending on machine)

c    2^31=2,147,483,648
      npart=5000000

      verbose=.false.
!      verbose=.true.
      iout=6
      nval=0
      sumtemp=0.0d0
      diffm=0.0d0

      call ppidd_size(nprocs)
      call ppidd_rank(iprocs)
! ARCCA Test Flag
      arcca_test_flag=.false.
!      arcca_test_flag=.true.

      num_tasks=50
      if(arcca_test_flag) num_tasks=1000*nprocs !ARCCA
      nscale_task=1


      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,*) 'PPIDD shared counter tests:'
      call ppidd_flush6

#if defined(GA_MPI)

      if(iprocs.eq.0) write(iout,*)
     > 'For GA_MPI, ppidd_nxtval does nothing. Skip this test.'

#else

c... shared counter test with computing tasks

c... initialise the buffer space
      ALLOCATE( task_array(nprocs) )
      ALLOCATE( sum_array(nprocs) )
      ALLOCATE( tcpu_array1(nprocs) )
      ALLOCATE( twt_array1(nprocs) )
      ALLOCATE( twt_array2(nprocs) )
      do i=1,nprocs
        task_array(i)=0
        sum_array(i)=0.0d0
        tcpu_array1(i)=0.0d0
        twt_array1(i)=0.0d0
        twt_array2(i)=0.0d0
      end do

c...(1)  shared counter test with computing tasks
      if(iprocs.eq.0) write(iout,*)
     >  '(1) PPIDD shared counter test with computing tasks:'
      if(verbose) write(iout,30)
     >            'Nprocs=',nprocs,'  Total tasks=', num_tasks
c... initialise the NXTVAL Server
      call ppidd_nxtval(-nprocs,junk)
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      call ppidd_nxtval(nprocs,nval)
      wtimea=ppidd_wtime()
      do while (nval.lt.num_tasks)
        call compute_task(npart,nscale_task,sumtemp)
        sum_array(iprocs+1)=sum_array(iprocs+1)+sumtemp
        task_array(iprocs+1)=task_array(iprocs+1)+1
        wtimeb=ppidd_wtime()
c... time spent on computation and other overhead
        twt_array2(iprocs+1) =twt_array2(iprocs+1)+wtimeb-wtimea
        call ppidd_nxtval(nprocs,nval)
        wtimea=ppidd_wtime()
      enddo
      call get_current_times(cpu2,wtime2)
c... task time
      tcpu_array1(iprocs+1)=cpu2-cpu1
      twt_array1(iprocs+1)=wtime2-wtime1
      call ppidd_barrier()


c... release the shared counter number in NXTVAL Server and set it to zero
      call ppidd_nxtval(-nprocs,junk)

c... broadcast the completed tasks on current process
      do i=1,nprocs
        call ppidd_bcast(c_loc(task_array(i)),1,0_c_int,i-1)
        call ppidd_bcast(c_loc(sum_array(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(tcpu_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array2(i)), 1,1_c_int,i-1)
      end do

c... check whether right numbers are obtained
      do i=2,nprocs
        diffm=sum_array(i)/dble(task_array(i))
     >        - sum_array(i-1)/dble(task_array(i-1))
        if (verbose) write(iout,*) ' difference of average=',diffm
        if (abs(diffm).gt.1.0d-5)
     >       call ppidd_error('shared counter test failed'//c_null_char,
     >                        0_c_int)
      end do

c... print out the summary
      if(iprocs.eq.0) then
        write(iout,30) 'Nprocs=',nprocs,'  Total tasks=', num_tasks
30      format(1x,a,t15,i4,a,t40,i12)
        write(iout,19)
19      format(1x,' IPROC      TASKS      WTime/task(Sec)')
        do i=1,nprocs
          if (task_array(i).ne.0) then
            twt_task=twt_array1(i)/dble(task_array(i))
          else
            twt_task=100000.0
          end if
          write(iout,20) i-1,task_array(i),twt_task
20        format(1x,i5,i12,f20.9)
        end do
        write(iout,10) 'Time spent for iprocs=0:',
     >  tcpu_array1(1),twt_array1(1)
10      format(1x,a,t49,' cpu=',f16.6,' sec,',
     >  '  wall time=',f16.6,' sec')
        cpu_tot1=sum(tcpu_array1)
        wt_tot1=sum(twt_array1)
        wt_tot2=sum(twt_array2)
        write(iout,10) 'Sum of time for all procs all tasks:',
     >  cpu_tot1,wt_tot1
        write(iout,9) 'Sum of time for all procs all computations:',
     >  wt_tot2
9       format(1x,a,t75,'  wall time=',f16.6,' sec')
        write(iout,9) 'Sum of time for all procs all counter callings:',
     >  wt_tot1-wt_tot2
        acpu1=cpu_tot1/dble(num_tasks)
        awt1=wt_tot1/dble(num_tasks)
        awt2=wt_tot2/dble(num_tasks)
        write(iout,10) 'Average time per task:',acpu1,awt1
        write(iout,9) 'Average time per computation:',awt2
        write(iout,12) 'Average time per shared counter calling:',
     >                 (awt1-awt2)*1.0d9
12      format(1x,a,t75,'  wall time=',f12.1,' nanosec')
      end if
c... end of (1)

c...(2)  distribute tasks evenly to every process without shared counter
      call ppidd_barrier()
      do i=1,nprocs
        task_array(i)=0
        sum_array(i)=0.0d0
      end do
      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,*) '(2) Distribute computing tasks to',
     >' processes without shared counter:'
      if(verbose) write(iout,30)
     >            'Nprocs=',nprocs,'  Total tasks=', num_tasks
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      do i=1,num_tasks
        ipr=mod(i,nprocs)
        if(ipr.eq.iprocs) then
          call compute_task(npart,nscale_task,sumtemp)
          sum_array(iprocs+1)=sum_array(iprocs+1)+sumtemp
          task_array(iprocs+1)=task_array(iprocs+1)+1
        end if
      enddo
      call get_current_times(cpu2,wtime2)
      tcpu_array1(iprocs+1)=cpu2-cpu1
      twt_array1(iprocs+1)=wtime2-wtime1
      call ppidd_barrier()

c... broadcast the completed tasks on current process
      do i=1,nprocs
        call ppidd_bcast(c_loc(task_array(i)),1,0_c_int,i-1)
        call ppidd_bcast(c_loc(sum_array(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(tcpu_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array1(i)), 1,1_c_int,i-1)
      end do

c... check whether right numbers are obtained
      do i=2,nprocs
        diffm=sum_array(i)/dble(task_array(i))
     >        - sum_array(i-1)/dble(task_array(i-1))
        if (verbose) write(iout,*) ' difference of average=',diffm
        if (abs(diffm).gt.1.0d-5)
     >       call ppidd_error('shared counter test failed'//c_null_char,
     >                        0_c_int)
      end do

c... print out the summary
      if(iprocs.eq.0) then
        write(iout,30) 'Nprocs=',nprocs,'  Total tasks=', num_tasks
        write(iout,19)
        do i=1,nprocs
          if (task_array(i).ne.0) then
            twt_task=twt_array1(i)/dble(task_array(i))
          else
            twt_task=100000.0
          end if
          write(iout,20) i-1,task_array(i),twt_task
        end do
        write(iout,10) 'Time spent for iprocs=0:',
     >  tcpu_array1(1),twt_array1(1)
        cpu_tot3=sum(tcpu_array1)
        wt_tot3=sum(twt_array1)
        write(iout,10) 'Sum of time for all procs all computing tasks:',
     >  cpu_tot3,wt_tot3
        acpu3=cpu_tot3/dble(num_tasks)
        awt3=wt_tot3/dble(num_tasks)
        write(iout,10) 'Average time per computing task:',
     >                 acpu3,awt3
        write(iout,*)
        write(iout,115)
115     format(1x,'Comparing the times between (1) and (2), we also'
     >  ,' get the approximate average time per counter calling:')
        write(iout,12) ' ',(awt1-awt3)*1.0d9
        write(iout,13)
13      format(1x,'This data probably can not be used to evaluate'
     >  ,' the shared counter since some extra overhead is included.')
      end if
c... end of (2) shared counter test

c... (3) shared counter test without computing tasks
c    2^31=2,147,483,648
      num_tasks=1000000
      if(arcca_test_flag) num_tasks=100000000/nprocs !ARCCA
      do i=1,nprocs
        task_array(i)=0
      end do
      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,*)
     >  '(3) PPIDD shared counter test without computing tasks:'
      if(verbose) write(iout,30)
     >            'Nprocs=',nprocs,'  Total tasks=', num_tasks
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      call ppidd_nxtval(nprocs,nval)
      do while (nval.lt.num_tasks)
        task_array(iprocs+1)=task_array(iprocs+1)+1
        call ppidd_nxtval(nprocs,nval)
      enddo
      call get_current_times(cpu2,wtime2)
      tcpu_array1(iprocs+1)=cpu2-cpu1
      twt_array1(iprocs+1)=wtime2-wtime1
      call ppidd_barrier()

c... release the shared counter number in NXTVAL Server and set it to zero
      call ppidd_nxtval(-nprocs,junk)

c... broadcast the completed tasks on current process
      do i=1,nprocs
        call ppidd_bcast(c_loc(task_array(i)),1,0_c_int,i-1)
        call ppidd_bcast(c_loc(tcpu_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array1(i)), 1,1_c_int,i-1)
      end do
      call ppidd_barrier()

c... print out the summary
      if(iprocs.eq.0) then
        write(iout,30) 'Nprocs=',nprocs,'  Total tasks=', num_tasks
        write(iout,19)
        do i=1,nprocs
          if (task_array(i).ne.0) then
            twt_task=twt_array1(i)/dble(task_array(i))
          else
            twt_task=100000.0
          end if
          write(iout,20) i-1,task_array(i),twt_task
        end do
        write(iout,10) 'Time spent for iprocs=0:',
     >  tcpu_array1(1),twt_array1(1)
        cpu_tot4=sum(tcpu_array1)
        wt_tot4=sum(twt_array1)
        write(iout,10) 'Sum of time for all procs all empty tasks:',
     >  cpu_tot4,wt_tot4
c... convert time to nanosecond
        acpu4=(cpu_tot4/dble(num_tasks))*1.0d9
        awt4=(wt_tot4/dble(num_tasks))*1.0d9
        write(iout,11) 'Average time per shared counter calling:',
     >                 acpu4,awt4
11      format(1x,a,t49,' cpu=',f12.1,' nanosec,',
     >  '  wall time=',f12.1,' nanosec')
        write(iout,111)
111     format(1x,'The data in (3) probably can not be used to evaluate'
     >  ,' the shared counter since all callings may concentrate on',
     >  ' some processes.')
        write(iout,*) 'The data in (1) are much better.'
      end if
c... end of (3) shared counter test

c... release the memory in buffer space
      IF( ALLOCATED(task_array) ) DEALLOCATE( task_array)
      IF( ALLOCATED(sum_array) ) DEALLOCATE( sum_array)
      IF( ALLOCATED(tcpu_array1) ) DEALLOCATE( tcpu_array1)
      IF( ALLOCATED(twt_array1) ) DEALLOCATE( twt_array1)
      IF( ALLOCATED(twt_array2) ) DEALLOCATE( twt_array2)

#endif

      return
      end



      subroutine ppidd_mutex_test(helper_server_flag)
      use ppidd
      implicit double precision (a-h,o-z)
      logical helper_server_flag
      integer,          allocatable, target :: task_array(:)
      double precision, allocatable, target :: sum_array(:)
      double precision, allocatable, target :: tcpu_array1(:)
      double precision, allocatable, target :: twt_array1(:)
      double precision cpu1,cpu2
      double precision wtime1,wtime2
      double precision sumtemp,diffm
      integer(c_int64_t) :: nprocs,iprocs
      integer iout,i
      integer num_tasks,nscale_task,npart
      integer(c_int64_t) :: storetype
      integer(c_int64_t) :: number
      integer(c_int)     :: ok
      logical verbose
      logical arcca_test_flag
c
c... set task number and task scale
c    please change these two numbers if carrying out large-scale, many-task test
c    eg. num_tasks=1000 nscale_task=1000 (Please be aware that it may take hours depending on machine)

c    2^31=2,147,483,648
      npart=5000000

      verbose=.false.
!      verbose=.true.
      iout=6
      nval=0
      sumtemp=0.0d0
      diffm=0.0d0

      call ppidd_size(nprocs)
      call ppidd_rank(iprocs)
! ARCCA Test Flag
      arcca_test_flag=.false.
!      arcca_test_flag=.true.
      num_tasks=50
      if(arcca_test_flag) num_tasks=1000 !ARCCA
      nscale_task=1

      num_tasks=(num_tasks/nprocs)*nprocs
      if(num_tasks.eq.0) num_tasks=nprocs

      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,*) 'PPIDD mutex tests:'
      call ppidd_flush6

c... calculate the average time per get_current_times calling
      ntimes=1000000
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      do i=1,ntimes
        call get_current_times(cpu0,wtime0)
      end do
      call get_current_times(cpu2,wtime2)
      cpu_tot=cpu2-cpu1
      wt_tot=wtime2-wtime1

      if(iprocs.eq.0) then
        write(iout,*)
        write(iout,10) 'Total time on calling get_current_times: ',
     >    cpu_tot,wt_tot
        write(iout,89) 'Total number of calling get_current_times: ',
     >    ntimes
89      format(1x,a,t49,i10)
c... convert time to nanosecond
        acpu=(cpu_tot/dble(ntimes))*1.0d9
        awt=(wt_tot/dble(ntimes))*1.0d9
        write(iout,11) 'Average time per get_current_times calling:',
     >     acpu,awt
      end if

c... calculate the average time per ppidd_wtime calling
      call ppidd_barrier()
      wtime1=ppidd_wtime()
      do i=1,ntimes
        wtime0=ppidd_wtime()
      end do
      wtime2=ppidd_wtime()
      wt_tot=wtime2-wtime1

      if(iprocs.eq.0) then
        write(iout,*)
        write(iout,92) 'Total time on calling ppidd_wtime ',
     >     ntimes, wt_tot
92      format(1x,a,i10,'  times:',t75,'  wall time=',f16.6,' sec')
c... convert time to nanosecond
        awt=(wt_tot/dble(ntimes))*1.0d9
        write(iout,12) 'Average time per ppidd_wtime calling:', awt
        write(iout,*)
      end if


c... initialise the buffer space
      ALLOCATE( task_array(nprocs) )
      ALLOCATE( sum_array(nprocs) )
      ALLOCATE( tcpu_array1(nprocs) )
      ALLOCATE( twt_array1(nprocs) )


c... create mutexes
      number=1
#ifdef MPI2
      if (helper_server_flag) then
        nloop_array=2
      else
        nloop_array=1
      end if
#else
      nloop_array=1
#endif
      do loop_array=1,nloop_array
      if (nloop_array.eq.1) then
        if(iprocs.eq.0) write(iout,101) ' '
101     format(/,a,'Mutex data are stored on compute process, and'
     >  ,' mutex is implemented through RMA operations:')
        storetype=0_c_int64_t
      else
c... case 1
        if (loop_array.eq.1) then
          if(iprocs.eq.0) write(iout,101) ' 1. '
          storetype=0_c_int64_t
c... case 2
        else
          if(iprocs.eq.0) write(iout,102)
102       format(/' 2. Mutex data are stored on helper process, and'
     >  ,' mutex is implemented on top of two-sided operations:')
          storetype=1_c_int64_t
        end if
      end if

      ok=ppidd_create_mutexes(storetype,number)

      do i=1,nprocs
        task_array(i)=0
        sum_array(i)=0.0d0
        tcpu_array1(i)=0.0d0
        twt_array1(i)=0.0d0
      end do

c...(1)  mutex test with computing tasks
      if(iprocs.eq.0) write(iout,*)
     >  '(1) PPIDD mutex test with computing tasks:'
      if(verbose) write(iout,30)
     >            'Nprocs=',nprocs,'  Total tasks=', num_tasks

      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)

      do i=iprocs+1,num_tasks,nprocs
        call ppidd_lock_mutex(0_c_int64_t)
c... critical section
        call compute_task(npart,nscale_task,sumtemp)
        sum_array(iprocs+1)=sum_array(iprocs+1)+sumtemp
        task_array(iprocs+1)=task_array(iprocs+1)+1
        call ppidd_unlock_mutex(0_c_int64_t)
        call ppidd_barrier()
      enddo

      call get_current_times(cpu2,wtime2)

c... task time
      tcpu_array1(iprocs+1)=cpu2-cpu1
      twt_array1(iprocs+1)=wtime2-wtime1
      call ppidd_barrier()

c... broadcast the completed tasks on current process
      do i=1,nprocs
        call ppidd_bcast(c_loc(task_array(i)),1,0_c_int,i-1)
        call ppidd_bcast(c_loc(sum_array(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(tcpu_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array1(i)), 1,1_c_int,i-1)
      end do

c... check whether right numbers are obtained
      do i=2,nprocs
        diffm=sum_array(i)/dble(task_array(i))
     >        - sum_array(i-1)/dble(task_array(i-1))
        if (verbose) write(iout,*) ' difference of average=',diffm
        if (abs(diffm).gt.1.0d-5)
     >       call ppidd_error('mutex test failed'//c_null_char,
     >                        0_c_int)
      end do

c... print out the summary
      if(iprocs.eq.0) then
        write(iout,30) 'Nprocs=',nprocs,'  Total tasks=', num_tasks
30      format(1x,a,t15,i4,a,t40,i12)
        write(iout,19)
19      format(1x,' IPROC      TASKS      WTime/task(Sec)')
        do i=1,nprocs
          if (task_array(i).ne.0) then
            twt_task=twt_array1(i)/dble(task_array(i))
          else
            twt_task=100000.0
          end if
          write(iout,20) i-1,task_array(i),twt_task
20        format(1x,i5,i12,f20.9)
        end do
        write(iout,10) 'Time spent for iprocs=0:',
     >  tcpu_array1(1),twt_array1(1)
10      format(1x,a,t49,' cpu=',f16.6,' sec,',
     >         '  wall time=',f16.6,' sec')
        cpu_tot1=sum(tcpu_array1)
        wt_tot1=sum(twt_array1)
        acpu_tot1=cpu_tot1/dble(nprocs)
        awt_tot1=wt_tot1/dble(nprocs)
        write(iout,10) 'Time spent for each proc(average):',
     >  acpu_tot1,awt_tot1
        acpu1=acpu_tot1/dble(num_tasks)
        awt1=awt_tot1/dble(num_tasks)
        write(iout,10) 'Average time per task:',acpu1,awt1
12      format(1x,a,t75,'  wall time=',f12.1,' nanosec')
      end if
c... end of (1)

c...(2)  distribute tasks evenly to every process without mutex (lock and unlock)
      call ppidd_barrier()
      do i=1,nprocs
        task_array(i)=0
        sum_array(i)=0.0d0
      end do
      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,23)
23    format(1x,'(2) Distribute computing tasks to',
     >' all processes without mutex, but ensure only one process is',
     >' allowed to do the computation at a time:')
      if(verbose) write(iout,30)
     >            'Nprocs=',nprocs,'  Total tasks=', num_tasks
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      do i=1,num_tasks
        if(mod(i,nprocs).eq.iprocs) then
          call compute_task(npart,nscale_task,sumtemp)
          sum_array(iprocs+1)=sum_array(iprocs+1)+sumtemp
          task_array(iprocs+1)=task_array(iprocs+1)+1
        end if
        call ppidd_barrier()
      enddo
      call get_current_times(cpu2,wtime2)
      tcpu_array1(iprocs+1)=cpu2-cpu1
      twt_array1(iprocs+1)=wtime2-wtime1
      call ppidd_barrier()

c... broadcast the completed tasks on current process
      do i=1,nprocs
        call ppidd_bcast(c_loc(task_array(i)),1,0_c_int,i-1)
        call ppidd_bcast(c_loc(sum_array(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(tcpu_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array1(i)), 1,1_c_int,i-1)
      end do

c... check whether right numbers are obtained
      do i=2,nprocs
        diffm=sum_array(i)/dble(task_array(i))
     >        - sum_array(i-1)/dble(task_array(i-1))
        if (verbose) write(iout,*) ' difference of average=',diffm
        if (abs(diffm).gt.1.0d-5)
     >       call ppidd_error('mutex test failed'//c_null_char,
     >                        0_c_int)
      end do

c... print out the summary
      if(iprocs.eq.0) then
        write(iout,30) 'Nprocs=',nprocs,'  Total tasks=', num_tasks
        write(iout,19)
        do i=1,nprocs
          if (task_array(i).ne.0) then
            twt_task=twt_array1(i)/dble(task_array(i))
          else
            twt_task=100000.0
          end if
          write(iout,20) i-1,task_array(i),twt_task
        end do
        write(iout,10) 'Time spent for iprocs=0:',
     >  tcpu_array1(1),twt_array1(1)
        cpu_tot3=sum(tcpu_array1)
        wt_tot3=sum(twt_array1)
        acpu_tot3=cpu_tot3/dble(nprocs)
        awt_tot3=wt_tot3/dble(nprocs)
        write(iout,10) 'Time spent for each proc(average):',
     >  acpu_tot3,awt_tot3
        acpu3=acpu_tot3/dble(num_tasks)
        awt3=awt_tot3/dble(num_tasks)
        write(iout,10) 'Average time per computing task:',
     >                 acpu3,awt3
        write(iout,*)
        write(iout,33)
33      format(1x, 'Comparing the time between (1) and (2), we also',
     >  ' get the approximate average time per lock/unlock calling:')
        write(iout,12) ' ',(awt1-awt3)*1.0d9
        write(iout,13)
13      format(1x,'This data probably can be used to evaluate'
     >  ,' the mutex only if the number of tasks is large enough.')
      end if
c... end of (2) mutex test

c... (3) mutex test without computing tasks
c    2^31=2,147,483,648
      num_tasks3=10000
      if(arcca_test_flag) num_tasks3=200000/nprocs !ARCCA
      do i=1,nprocs
        task_array(i)=0
      end do
      if(iprocs.eq.0) write(iout,*)
      if(iprocs.eq.0) write(iout,*)
     >  '(3) PPIDD mutex test without computing tasks:'
      if(verbose) write(iout,30)
     >            'Nprocs=',nprocs,'  Total tasks=', num_tasks3
      call ppidd_barrier()
      call get_current_times(cpu1,wtime1)
      do i=1,num_tasks3
        call ppidd_lock_mutex(0_c_int64_t)
        call ppidd_unlock_mutex(0_c_int64_t)
      enddo
      call get_current_times(cpu2,wtime2)
      call ppidd_barrier()
      call get_current_times(cpu3,wtime3)
      task_array(iprocs+1)=num_tasks3
      tcpu_array1(iprocs+1)=cpu2-cpu1
      twt_array1(iprocs+1)=wtime2-wtime1


c... broadcast the completed tasks on current process
      do i=1,nprocs
        call ppidd_bcast(c_loc(task_array(i)),1,0_c_int,i-1)
        call ppidd_bcast(c_loc(tcpu_array1(i)), 1,1_c_int,i-1)
        call ppidd_bcast(c_loc(twt_array1(i)), 1,1_c_int,i-1)
      end do
      call ppidd_barrier()

c... print out the summary
      num_tasks_tot=sum(task_array)
      if(iprocs.eq.0) then
        write(iout,30) 'Nprocs=',nprocs,'  Total tasks=', num_tasks_tot
        write(iout,19)
        do i=1,nprocs
          if (task_array(i).ne.0) then
            twt_task=twt_array1(i)/dble(task_array(i))
          else
            twt_task=100000.0
          end if
          write(iout,20) i-1,task_array(i),twt_task
        end do
        write(iout,10) 'Time spent for iprocs=0:',
     >  tcpu_array1(1),twt_array1(1)
        cpu_tot4=sum(tcpu_array1)
        wt_tot4=sum(twt_array1)
        acpu_tot4=cpu_tot4/dble(nprocs)
        awt_tot4=wt_tot4/dble(nprocs)
        write(iout,10) 'Time spent for each proc(average):',
     >  acpu_tot4,awt_tot4
c... convert time to nanosecond
        acpu4=(acpu_tot4/dble(num_tasks_tot))*1.0d9
        awt4=(awt_tot4/dble(num_tasks_tot))*1.0d9

        write(iout,11) 'Average time per lock/unlock calling:',
     >                 acpu4,awt4
11      format(1x,a,t49,' cpu=',f12.1,' nanosec,',
     >  '  wall time=',f12.1,' nanosec')
        write(iout,111)
111     format(1x,'The data in (3) can be used to evaluate the mutex'
     >  ,' since all callings are evenly distributed to all processes.')
!        write(iout,*) 'The data in (1) are much better.'
      end if
c... end of (3) mutex test

c... destroy mutexes
      ok=ppidd_destroy_mutexes()

      end do
c...  end of loop_array

c... release the memory in buffer space
      IF( ALLOCATED(task_array) ) DEALLOCATE( task_array)
      IF( ALLOCATED(sum_array) ) DEALLOCATE( sum_array)
      IF( ALLOCATED(tcpu_array1) ) DEALLOCATE( tcpu_array1)
      IF( ALLOCATED(twt_array1) ) DEALLOCATE( twt_array1)

      return
      end



      subroutine compute_task(N,nloop,pai)
c
c This simple subroutine approximates pi by computing pi = integral
c from 0 to 1 of 4/(1+x*x)dx which is approximated by sum from
c k=1 to N of 4 / ((1 + (k-1/2)**2 ).  The input data required are N and nloop.
c nloop is used only for duplicating computation
c The returned value sum is computed pi.
c
      implicit double precision (a-h,o-z)
      integer N,nloop
      double precision totsum,sum,pai
      integer i,j
      double precision err, pi, w
      pi = 4.0d0*atan(1.0d0)

      totsum = 0.0d0
      err = 0.0d0
      do i=1,nloop
        sum = 0.0d0
        if (N .gt. 0) then
          w = 1.0d0/dble(N)
          do j = 1,N
             sum = sum + f((dble(j)-0.5d0)*w)
          enddo
          sum = sum * w
        else
          write(6,*) 'ERROR in running compute_task: N .le. 0. STOP'
          stop
        end if
        totsum=totsum+sum
      end do
      pai=totsum/dble(nloop)
      err = pai - pi
!      write(6,*) 'running compute_task, pai=',pai,' err=',err

      return
      contains
       double precision function f(x)
       implicit none
       double precision x
       f = 4.0d0/(1.0d0+x*x)
       return
       end function f

      end
