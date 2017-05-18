   program mpactread
   use  hdf5
   use  mod_hdftools
   use  mod_batch
   implicit none
!=======================================================================
!
!  Program to read MPACT HDF output file and print summary
!
!  Copyright (c) 2014-2017 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2014/04/28 - draft 1
!  2014/05/01 - update for VERAOUT version 1
!  2014/08/01 - change HDF open statement to read-only
!             - read exposure values
!  2015/03/06 - add additional scalar edits
!             - added option to print single distributions with "d" command line arguments
!  2015/03/13 - set fuel temperatures in guide tubes to zero
!             - calculate averages in qtr-symmetry problems correctly
!             - add pin exposure edits
!  2015/03/25 - add outer iterations count and time edits
!  2017/04/25 - add initial batch edit option for pin exposures
!             - assumes assm_map is 2 character labels ***
!  2017/05/18 - add 2PIN and 2PIN map for Jim
!
!-----------------------------------------------------------------------

      character(len=80)  :: inputfile
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i, j, k, n
      integer            :: idis
      integer            :: ierror
      integer            :: itype
      integer            :: ndim         ! temp variable
      integer            :: idim(10)     ! temp variable
      integer            :: nstate=0     ! statepoint number

      integer            :: llpow       ! index for power array
      integer            :: llexp       ! index for exposure array
      integer            :: lltfu       ! index for fuel temperature

      logical            :: ifxst
      logical            :: ifdebug=.false. ! debug flag
      logical            :: ifoldstate=.false.

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=12)  :: group_name

      integer(hid_t)     :: file_id

! input data

      integer  :: isym                      ! core symmetry
      integer  :: iver                      ! output file version number
      integer  :: naxial                    ! number of axial levels
      integer  :: icore, jcore              ! core size
      integer, allocatable :: mapcore(:,:)  ! core map

      integer  :: kd                        ! number of axial levels in pin maps
      integer  :: nassm                     ! number of assemblies in pin maps
      integer  :: npin                      ! pin data

      real(8)  :: xtemp                     ! temp value
      real(8)  :: xtemp2                    ! temp value
      real(8)  :: xave                      ! temp value
      real(8)  :: xkeff                     ! eigenvalue
      real(8)  :: xexpo                     ! exposure
      real(8)  :: xefpd                     ! exposure EFPD
      real(8)  :: xflow                     !
      real(8)  :: xpow                      !
      real(8)  :: boron                     !
      real(8)  :: rated_power               ! rated flow
      real(8)  :: rated_flow                ! rated power
      real(8)  :: apitch                    ! assembly pitch
      real(8)  :: outlet_temp               ! outlet temperature

      real(8), allocatable :: axial(:)      ! axial elevations
      real(8), allocatable :: temp4d(:,:,:,:)
      real(8), allocatable :: power(:,:,:,:)  ! 3d distribution
      real(8), allocatable :: tdist(:,:,:,:)  ! 3d distribution
      real(8), allocatable :: tdist2d(:,:,:)  ! 2d collapsed distribution

      character(len=80)  :: title           ! Problem title
      character(len=2), allocatable :: xlabel(:)  ! assembly map labels
      character(len=2), allocatable :: ylabel(:)  ! assembly map labels

      integer, parameter :: maxdist=6    ! maximum number of distributions
      character(len=20) :: dist_label(maxdist)
      logical           :: dist_print(maxdist)   ! logical to print distribution

      integer, parameter :: maxstate=200
      real(8)  :: state_xkeff(maxstate)
      real(8)  :: state_xexpo(maxstate)  ! exposure
      real(8)  :: state_xefpd(maxstate)
      real(8)  :: state_boron(maxstate)
      real(8)  :: state_flow (maxstate)
      real(8)  :: state_power(maxstate)
      real(8)  :: state_tinlet(maxstate)
      real(8)  :: state_3exp(maxstate)
      real(8)  :: state_3pin(maxstate)
      real(8)  :: state_2pin(maxstate)
      integer  :: state_nout(maxstate)
      real(8)  :: state_time(maxstate)

! command line flags

      logical  :: if3d   =.false.           ! turn on 3D edits
      logical  :: if2d   =.false.           ! turn on 2D edits
      logical  :: if2da  =.false.           ! turn on 2D assembly edits
      logical  :: if2pin =.false.           ! turn on 2PIN assembly edits
      logical  :: if1d   =.false.           ! turn on 1D edits
      logical  :: iftime =.false.           ! timing summary
      logical  :: ifbatch=.false.           ! batch edits

!  initialize

      inputfile=' '

      state_xkeff(:)=0.0d0
      state_xexpo(:)=0.0d0
      state_xefpd(:)=0.0d0
      state_boron(:)=0.0d0
      state_flow (:)=0.0d0
      state_power(:)=0.0d0
      state_tinlet(:)=-1.0d0
      state_3exp(:)=0.0d0
      state_3pin(:)=0.0d0
      state_2pin(:)=0.0d0
      state_nout(:)=0
      state_time(:)=-1.0d0

      dist_label(1)='pin_powers'     ! must come first
      dist_label(2)='pin_fueltemps'
      dist_label(3)='pin_cladtemps'
      dist_label(4)='pin_modtemps'
      dist_label(5)='pin_moddens'
      dist_label(6)='pin_exposures'

      llpow=1       ! save location of power
      lltfu=2       ! save location of fuel temperatures
      llexp=6       ! save location of exposures

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  mpactread.exe [hdf5_file] {1D/2D/2DA/3D} {time} {-dN}'
        write (*,*) '  1D     print 1D edits'
        write (*,*) '  2D     print 2D pin edits'
        write (*,*) '  2DA    print 2D assembly edits'
        write (*,*) '  3D     print 3D pin edits'
        write (*,*) '  time   print timing summary'
        write (*,*) '  batch  print batch edits'
        write (*,*)
        write (*,*) 'distribution selections for -dN option:'
        do idis=1, maxdist
          write (*,18) idis, trim(dist_label(idis))
        enddo
        stop
      endif
   18 format ('   -d',i1,5x,a)

! parse command line arguments

      idis=-1
      dist_print(:)=.false.

      do i=1, iargs
        call get_command_argument(i,carg)
        if     (carg.eq.'3D' .or. carg.eq.'3d') then
          if3d=.true.
        elseif (carg.eq.'2D' .or. carg.eq.'2d') then
          if2d=.true.
        elseif (carg.eq.'2DA' .or. carg.eq.'2da') then
          if2da=.true.
        elseif (carg.eq.'2PIN' .or. carg.eq.'2pin') then
          if2pin=.true.
        elseif (carg.eq.'1D' .or. carg.eq.'1d') then
          if1d=.true.
        elseif (carg.eq.'debug') then
          ifdebug=.true.
        elseif (carg.eq.'time') then
          iftime=.true.
        elseif (carg.eq.'batch') then
          ifbatch=.true.
        elseif (carg.eq.'-help' .or. carg.eq.'--help') then  ! add for Ben
          write (*,*) 'run mpactread with no command line arguments for help'
        elseif (carg(1:2).eq.'-d') then
          read (carg(3:),*) idis
          if (idis.ge.1 .and. idis.le.maxdist) dist_print(idis)=.true.
        else
          inputfile=carg
        endif
      enddo

      if (idis.eq.-1) then   ! print all distributions by default
        dist_print(:)=.true.
      endif

      if (if1d)   write (*,*) '1D edits requested on command line'
      if (if2d)   write (*,*) '2D pin edits requested on command line'
      if (if2da)  write (*,*) '2DA assembly edits requested on command line'
      if (if2pin) write (*,*) '2PIN edits requested on command line'
      if (if3d)   write (*,*) '3D pin edits requested on command line'

      if (inputfile.eq.' ') then
        stop 'no input file specified on command line'
      endif

!--- initialize fortran interface

      call h5open_f(ierror)

!--- open HDF file

      write (*,'(2a)') 'reading h5 file: ', trim(inputfile)
      inquire(file=inputfile, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(inputfile),' does not exist'
        stop
      endif
      call h5fopen_f (inputfile, H5F_ACC_RDONLY_F, file_id, ierror)   ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop
      endif

!--------------------------------------------------------------------------------
! Read top level HDF data
!--------------------------------------------------------------------------------

      xkeff=-100.0d0
      xexpo=-100.0d0
      iver=-100

! version

      dataset='veraout_version'
      call hdf5_read_integer(file_id, dataset, iver)

! title

      title=' '

      dataset='title'
      call read_string(file_id, dataset, title)

      write (*,*)
      write (*,'(2a)')     ' Title: ', trim(title)
      write (*,'(a,i2)') ' VERAOUT version ', iver
      write (*,*)

!-------------------
!  Read CORE group
!-------------------

      isym=-100
      idim(:)=0
      naxial=0          ! number of axial edit bounds
      icore=0           ! number of assemblies across
      jcore=0           ! number of assemblies across

      apitch=21.5d0     ! default value
      rated_power=0.0d0 ! default value
      rated_flow =0.0d0 ! default value

! check of CORE group exists (group may not exist on old MPACT files and all data will be in root)

      group_name='/CORE/'
      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (.not.ifxst) then
        group_name=' '
        write (*,*) 'CORE group not found - is this an old MPACT file?'
      endif

!--- assembly pitch

      dataset=trim(group_name)//'apitch'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, apitch)
      else
        write (*,*) '**** assembly pitch is missing from file ****'
      endif

!--- symmetry

      dataset=trim(group_name)//'core_sym'
      call hdf5_read_integer(file_id, dataset, isym)

!--- rated power and flow

      dataset=trim(group_name)//'rated_power'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, rated_power)
      else
        write (*,*) '**** rated power is missing from file ****'
      endif

      dataset=trim(group_name)//'rated_flow'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, rated_flow)
      else
        write (*,*) '**** rated flow  is missing from file ****'
      endif

      write (*,*)
      write (*,'(a,i2)') ' core symmetry  isym  = ', isym
      write (*,30)       'rated power    ', rated_power,' MW'
      write (*,30)       'rated flow     ', rated_flow
      write (*,30)       'assembly pitch ', apitch,' cm'
  30  format (1x,a,f12.4,a)

!--- core map

      dataset=trim(group_name)//'core_map'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.2) stop 'invalid dimensions in core_map'
      icore=idim(1)
      jcore=idim(2)
      if (ifdebug) write (*,*) 'debug: icore = ', icore
      if (ifdebug) write (*,*) 'debug: jcore = ', jcore
      allocate (mapcore(icore,jcore))
      call hdf5_read_integer(file_id, dataset, icore, jcore, mapcore)

      write (*,'(/,a)') ' Core Map:'
      do j=1, jcore
         write (*,'(2x,20i3)') (mapcore(i,j),i=1,icore)
      enddo

!--- axial mesh

      dataset=trim(group_name)//'axial_mesh'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.1) stop 'invalid dimensions in axial_mesh'
      naxial=idim(1)   ! actually naxial+1
      if (ifdebug) write (*,*) 'debug: naxial+1 = ', naxial
      allocate (axial(naxial))
      call hdf5_read_double(file_id, dataset, naxial, n, axial)   ! n will be set to naxial
      if (n.ne.naxial) then
         write (*,*) 'found    ', n
         write (*,*) 'expecting', naxial+1
         stop 'error reading axial'
      endif
      naxial=naxial-1    ! decrease to match number of levels instead of boundaries

      if (naxial.eq.0) then   ! fixup OLD MPACT
        write (*,*) 'ERROR in axial edits - is this an old broken MPACT file?'
        naxial=1
        axial(1)=0.0d0
        axial(2)=1.0d0
      endif

      write (*,'(/,a)') ' Axial Edit Boundaries'
      do j=1, naxial+1
        write (*,'(i5,20f12.4)') j-1, axial(j)
      enddo
      if (naxial.eq.0) then
        write (*,*) '*** warning: no axial edits found ****'
      endif

      write (*,'(a,f12.4)') ' total core height = ', axial(naxial+1)-axial(1)
      write (*,*)

      do j=1, naxial
        axial(j)=axial(j+1)-axial(j)   ! convert to deltas
      enddo
      axial(naxial+1)=0.0d0

      if (ifdebug) then
        write (*,*) 'axial deltas'
        do j=1, naxial+1
          write (*,'(i5,20f12.4)') j, axial(j)
        enddo
      endif

!--- temp define edit labels

  !  xlabel  R P N M L K J H G  F  E  D  C  B  A
  !  ylabel  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

      allocate (xlabel(icore))
      allocate (ylabel(jcore))

!**** still having issues with reading strings

!     dataset=trim(group_name)//'xlabel'
!     call h5info(file_id, dataset, itype, ndim, idim)
!     if (ndim.ne.1)        stop 'invalid dimensions in xlabel'
!     if (idim(1).ne.jcore) stop 'invalid size in xlabel'
!     call read_string1d(file_id, dataset, xlabel, icore)

!     dataset=trim(group_name)//'ylabel'
!     call h5info(file_id, dataset, itype, ndim, idim)
!     if (ndim.ne.1)        stop 'invalid dimensions in ylabel'
!     if (idim(1).ne.jcore) stop 'invalid size in ylabel'
!     call read_string1d(file_id, dataset, ylabel, jcore)

      n=ichar('A')+icore-1
      if (icore.ge.15) n=n+1   ! skip letter "O"
      if (icore.ge.9)  n=n+1   ! skip letter "I"
      do i=1, icore
        if (achar(n).eq.'I') n=n-1
        if (achar(n).eq.'O') n=n-1
        write (xlabel(i),'(a1,1x)') achar(n)
        n=n-1
      enddo
      do j=1, jcore
        write (ylabel(j),'(i2.2)') j
      enddo

!--- initialize batch edits

      if (ifbatch) then
        call batch_init(file_id, mapcore, icore, jcore)
      endif

!---------------------
!  Read STATE groups
!---------------------

!--- check if this file uses an old mpact STATE format

      group_name='/STATE_1/'
      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (ifxst) ifoldstate=.true.

! 40 format (/,'--------------------------',&
!            /,'  Reading statepoint ', i0, &
!            /,'--------------------------')

      nstate=0
      do
        nstate=nstate+1

        if (ifoldstate) then
          write (group_name,'(a,i0,a)') '/STATE_', nstate, '/'
        else
          write (group_name,'(a,i4.4,a)') '/STATE_', nstate, '/'
        endif
        if (ifdebug) write (*,*) 'debug: state= ', group_name

!--- check if statepoint exists

        call h5lexists_f(file_id, group_name, ifxst, ierror)
        if (.not.ifxst) then
          if (ifdebug) write (*,*) 'dataset not found - exiting'
          goto 800
        endif

!!      write (*,40) nstate   ! statepoint is already printed in reading dataset edits

        write (*,*)

!--------------------------------------------------------------------------------
! Read scalars
!--------------------------------------------------------------------------------

        xexpo=-1.0d0   ! initialize in case it is missing
        xefpd=-1.0d0
        xkeff=-1.0d0
        boron=-1.0d0
        xflow=-1.0d0
        xpow =-1.0d0

        dataset=trim(group_name)//'exposure'
        call hdf5_read_double(file_id, dataset, xexpo)

        dataset=trim(group_name)//'exposure_efpd'
        call hdf5_read_double(file_id, dataset, xefpd)

        dataset=trim(group_name)//'keff'
        call hdf5_read_double(file_id, dataset, xkeff)

        dataset=trim(group_name)//'boron'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, boron)
        endif

        dataset=trim(group_name)//'flow'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xflow)
        endif

        dataset=trim(group_name)//'power'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xpow)
        endif

        dataset=trim(group_name)//'tinlet'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, state_tinlet(nstate))
        endif

        dataset=trim(group_name)//'outer_timer'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, state_time(nstate))
        endif

        dataset=trim(group_name)//'outers'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_integer(file_id, dataset, state_nout(nstate))
        endif

        dataset=trim(group_name)//'outlet_temp'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, outlet_temp)
          write (*,*) 'outlet_temp = ', outlet_temp
        endif

        if (ifdebug) write (*,*) 'debug: keff = ', xkeff

!--- save statepoint values

        if (nstate.gt.maxstate) stop 'maxstate exceeded - increase and recompile'

        state_xkeff(nstate)=xkeff
        state_xexpo(nstate)=xexpo
        state_xefpd(nstate)=xefpd
        state_boron(nstate)=boron
        state_flow (nstate)=xflow
        state_power(nstate)=xpow

        write (*,'(a, f10.4,a)') ' exposure =', xexpo,' GWD/MT'
        write (*,'(a, f10.4,a)') ' exposure =', xefpd,' EFPD'
        write (*,'(a, f12.7)')   ' keff     =', xkeff

!--------------------------------------------------------------------------------
! Read pin powers
!--------------------------------------------------------------------------------

!--- allocate arrays if first statepoint

        if (nstate.eq.1) then
          dataset=trim(group_name)//'pin_powers'
          call h5info(file_id, dataset, itype, ndim, idim)

          nassm=idim(1)    ! order (nassm, kd, npin, npin)
          kd   =idim(2)
          npin =idim(4)
          if (ifdebug) then
            write (*,*) 'debug pin: ndim  = ', ndim
            write (*,*) 'debug pin: nassm = ', nassm
            write (*,*) 'debug pin: kd    = ', kd
            write (*,*) 'debug pin: npin  = ', npin
          endif

          if (ndim.ne.4)  stop 'invalid dimensions in pin data'
          if (nassm.eq.0) stop 'invalid number of assemblies in pin data'
          if (npin .eq.0) stop 'invalid number of pins in pin data'
          if (idim(3).ne.idim(4)) stop 'invalid npin in pin data'

          if (naxial.ne.kd) then
            stop 'mismatch found between axial levels in input and number of levels of output'
          endif

          if (ifdebug) write (*,*) 'debug: allocating pin arrays'

          allocate (temp4d(nassm, kd, npin, npin))
          allocate (power    (npin, npin, kd, nassm))    ! save power
          allocate (tdist    (npin, npin, kd, nassm))    ! 3D distributions
          allocate (tdist2d  (npin, npin, nassm))        ! 2D distributions
        endif

!-- Loop over distributions

        do idis=1, maxdist

          if (.not.dist_print(idis) .and. idis.ne.llpow) cycle    ! skip this distribution

          write (*,*)

          dataset=trim(group_name)//trim(dist_label(idis))
          call h5lexists_f(file_id, dataset, ifxst, ierror)
          if (.not.ifxst) then
            write (*,'(2a)') '>>', trim(dataset)
            write (*,'(2a)') 'WARNING: dataset does not exist for this file'
            cycle
          endif

          call hdf5_read_double(file_id, dataset, nassm, kd, npin, npin, temp4d)

          do n=1, nassm
            do k=1, kd
              do j=1, npin
                do i=1, npin
                  tdist(i,j,k,n)=temp4d(n,k,j,i)
                enddo
              enddo
            enddo
          enddo

          if (idis.eq.llpow) then
            power=tdist        ! save 3D power - needed for mask
          endif
          if (.not.dist_print(idis)) cycle    ! skip if this was power

          if (idis.eq.lltfu) then    ! apply mask to fuel temperatures
            call masktfu(npin,  kd, nassm, power, tdist)
          endif

          call stat3d(dist_label(idis), npin,  kd, nassm, icore, jcore, mapcore, &
                      axial, tdist, xave, xtemp, xtemp2)
          if (idis.eq.llpow) then
            state_3pin(nstate)=xtemp
            state_2pin(nstate)=xtemp2
          endif
          if (idis.eq.llexp) then
            state_3exp(nstate)=xtemp
          endif

          if (if2d .or. if2da) then
            call collapse2d(dist_label(idis), npin, kd, nassm, axial, tdist, tdist2d)
          endif

          if (idis.eq.llpow) then    ! special edits just for power
            call check_normalization(npin, kd, nassm, tdist, axial, icore, jcore, mapcore)
          endif

          if (if3d) then
            call print_3D_pin_map(dist_label(idis), npin, kd, nassm, tdist)
          endif
          if (if2d) then
            call print_3D_pin_map('2D '//dist_label(idis), npin, 1,  nassm, tdist2d)
          endif
          if (if1d) then
            call print1d(dist_label(idis), npin, kd, nassm, tdist, axial)
          endif

          if (if2da) then     ! 2D average assembly edits (already collapsed)
            call print2d_assm_map(dist_label(idis), npin, nassm, tdist2d, &
                 icore, jcore, mapcore, xlabel, ylabel)
          endif

          if (if2pin) then
            call print2d_2pin_map(dist_label(idis), npin, nassm, tdist, axial, &
                 icore, jcore, kd, mapcore, xlabel, ylabel)
          endif

          if (ifbatch .and. idis.eq.llexp) then
             call batchstat (npin, kd, nassm, icore, jcore, mapcore, axial, tdist, power)
          endif

        enddo   ! distributions

        if (ifbatch) then    ! batch edits
          call batchedit
        endif

      enddo   ! end of statepoint loop

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------
  800 continue

      call h5fclose_f(file_id, ierror)

!--- deallocate distributions

      if (allocated(tdist2d)) then  ! protect from missing statepoints
        deallocate (tdist2d)
        deallocate (tdist)
        deallocate (power)
        deallocate (temp4d)
      endif

      deallocate (mapcore)
      deallocate (xlabel,ylabel)
      deallocate (axial)

      nstate=nstate-1   ! decrease due to statepoint check



!--- print summary

      write (*,110)
      do n=1, nstate
        write (*,120) n, state_xexpo(n), state_xefpd(n), state_xkeff(n), &
               state_boron(n), state_3pin(n), state_2pin(n), state_3exp(n), state_flow(n), &
               state_power(n), state_tinlet(n)
      enddo
  110 format (/,'==================================',&
              /,'       Statepoint Summary', &
              /,'==================================',&
              /,'   N   exposure  exposure  eigenvalue   boron      3PIN      2PIN      3EXP      flow     power    tinlet')

  120 format (i4, f10.4, f10.2, f12.6, f10.2, 6f10.4)

!--- timing summary

      if (iftime) then
        write (*,150)
        do n=1, nstate
          write (*,160) n, state_xexpo(n), state_xefpd(n), state_nout(n), state_time(n)
        enddo
      endif
  150 format (/,'==================================',&
              /,'       Timing Summary', &
              /,'==================================',&
              /,'   N   exposure  exposure  outers   time (sec)')
  160 format (i4, f10.4, f10.2, i8, 1x, f12.3)

!--- finished

      write (*,'(/,a)') 'done'

      end program

!=======================================================================
!
!  Debug subroutine to check power normalization
!
!=======================================================================
      subroutine check_normalization(npin, kd, nassm, power, axial, icore, jcore, mapcore)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      integer, intent(in) :: icore, jcore
      integer, intent(in) :: mapcore(icore,jcore)
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)

!--- local

      integer :: i, j, k
      integer :: ia, ja, na
      real(8) :: pp
      real(8) :: cave, clen

      cave=0.0d0         ! core average
      clen=0.0d0         ! core fuel length

      write (*,'(/,1x,a)') 'Checking power normalization'

!--- check 3D normalization - loop over mapcore in case symmetry was used

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        do j=1, npin
          do i=1, npin
            do k=1, kd   ! loop over axial levels
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                cave=cave+axial(k)*pp
                clen=clen+axial(k)
              endif
            enddo
          enddo
        enddo

      enddo      ! ia
      enddo      ! ja

      if (clen.gt.0.0d0) cave=cave/clen

      if (abs(cave-1.0d0).gt.1.0d-6) then
         write (*,*) 'WARNING: pin powers do not appear to be normalized correctly'
         write (*,*) 'total core length ', clen
         write (*,*) 'core average power of linear power ', cave
      else
         write (*,*) 'core average  = ', cave
         write (*,*) 'normalization is correct'
      endif

      return
      end subroutine check_normalization
!=======================================================================
!
!  MPACT has non-zero values for fuel temperatures so zero the values
!  out if the power is zero
!
!=======================================================================
      subroutine masktfu (npin, kd, nassm, power, tdist)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8)             :: tdist(npin,npin,kd,nassm)

!--- local

      integer :: i, j, k, na
      real(8) :: pp

      write (*,'(/,1x,a)') 'Applying power mask to fuel temperature array'

!--- apply mask

      do na=1, nassm
        do k=1, kd   ! loop over axial levels
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.eq.0.0d0) then
                tdist(i,j,k,na)=0.0d0
              endif
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine masktfu

