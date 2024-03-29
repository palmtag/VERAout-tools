   program mpactread
   use  hdf5
   use  mod_hdftools
   use  mod_batch
   use  mod_coregeom, only : readcore, mapcore, icore, jcore, axial, nassm, npin, kd, pinload
   use  mod_nea,      only : nea_open, print_nea, nea_close
   implicit none
!=======================================================================
!
!  Program to read MPACT HDF output file and print summary
!
!  Copyright (c) 2014-2020 Core Physics, Inc.
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
!             - add ability to just print edits for a single statepoint
!             - add axial offset edit
!  2018/02/01 - added more control over pin loading edits (2d, 3d, 1d, etc.)
!  2018/04/26 - added csv option to print summary to csv file
!  2020/04/16 - update distribution names
!  2020/04/19 - add fuel temperature edits (-tfu)
!  2020/12/28 - add BWR T/H edits (-bwr)
!
!-----------------------------------------------------------------------

      character(len=80)  :: inputfile
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i, j, k, n
      integer            :: idis
      integer            :: ierror
      integer            :: idim(10)        ! temp variable
      integer            :: nstate=0        ! statepoint number
      integer            :: istate=0        ! user specified single statepoint option

      integer            :: llpow           ! index for power array
      integer            :: llexp           ! index for exposure array
      integer            :: lltfu           ! index for fuel temperature

      logical            :: ifxst
      logical            :: iflag
      logical            :: ifdebug=.false. ! debug flag
      logical            :: ifload =.false. ! print pin loadings

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=12)  :: group_name

      integer(hid_t)     :: file_id

! input data

      integer  :: iver                      ! output file version number

      real(8)  :: xtemp                     ! temp value
      real(8)  :: xtemp2                    ! temp value
      real(8)  :: xave                      ! temp value
      real(8)  :: xkeff                     ! eigenvalue
      real(8)  :: pressure                  ! pressure (MPa)

      real(8), allocatable :: temp4d(:,:,:,:)
      real(8), allocatable :: power(:,:,:,:)  ! 3d distribution
      real(8), allocatable :: tdist(:,:,:,:)  ! 3d distribution
      real(8), allocatable :: tdist2d(:,:,:)  ! 2d collapsed distribution

      character(len=80)  :: title              ! Problem title

      integer, parameter :: maxdist=10         ! maximum number of distributions
      character(len=30) :: dist_label(maxdist) ! distribution labels
      logical           :: dist_print(maxdist) ! logical to print distribution

! arrays for statepoint summary

      integer :: maxstate
      real(8), allocatable  :: state_xkeff(:)
      real(8), allocatable  :: state_expo(:)  ! exposure
      real(8), allocatable  :: state_efpd(:)  ! exposure EFPD
      real(8), allocatable  :: state_boron(:)
      real(8), allocatable  :: state_flow (:)
      real(8), allocatable  :: state_power(:)
      real(8), allocatable  :: state_tinlet(:)
      real(8), allocatable  :: state_3exp(:)
      real(8), allocatable  :: state_3pin(:)
      real(8), allocatable  :: state_2pin(:)
      real(8), allocatable  :: state_axoff(:)
      integer, allocatable  :: state_nout(:)
      real(8), allocatable  :: state_time(:)
      real(8), allocatable  :: state_atfu(:)
      real(8), allocatable  :: state_3tfu(:)

! command line flags

      logical  :: if3d   =.false.           ! turn on 3D edits
      logical  :: if2d   =.false.           ! turn on 2D edits
      logical  :: if2da  =.false.           ! turn on 2D assembly edits
      logical  :: if2pin =.false.           ! turn on 2PIN assembly edits
      logical  :: if1d   =.false.           ! turn on 1D edits
      logical  :: ifexit =.false.           ! turn on exit maps
      logical  :: iftime =.false.           ! timing summary
      logical  :: iftfu  =.false.           ! tfu edits
      logical  :: ifbwr  =.false.           ! bwr edits (t/h edits)
      logical  :: ifbatch=.false.           ! batch edits
      logical  :: ifcsv  =.false.           ! print summary to csv file
      logical  :: ifnea  =.false.           ! print summary to csv file in NEA benchmark format

!--- initialize

      inputfile=' '

      dist_label(1)='pin_powers'     ! pin power must come first in list
      dist_label(2)='pin_fuel_temp'         ! update 4/2020 VERA4.0
      dist_label(3)='pin_max_clad_surface_temp'     ! update 4/2020 VERA4.0
      dist_label(4)='pin_mod_temps'         ! update 4/2020 VERA4.0
      dist_label(5)='pin_mod_dens'          ! update 4/2020 VERA4.0
      dist_label(6)='pin_exposures'
      dist_label(7)='pin_steamrate'         ! add    4/2020 VERA4.0
      dist_label(8)='pin_cool_dens'         ! add 12/10/2020 BWR VERA4.2
      dist_label(9)='pin_cool_temp'         ! add 12/10/2020 BWR VERA4.2
      dist_label(10)='pin_cool_void'        ! add 12/10/2020 BWR VERA4.2

      llpow=1       ! save location of power
      lltfu=2       ! save location of fuel temperatures
      llexp=6       ! save location of exposures

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  mpactread.exe [hdf5_file] {options}'
        write (*,*) 'options:'
        write (*,*) '  1D     print 1D edits'
        write (*,*) '  2D     print 2D pin edits'
        write (*,*) '  2DA    print 2D assembly average edits'
        write (*,*) '  2PIN   print max 2D assembly rod edits'
        write (*,*) '  3D     print 3D pin edits'
        write (*,*) '  -exit  print exit maps for top plane'
        write (*,*) '  -load  print pin loadings'
        write (*,*) '  -csv   write summary to CSV file'
        write (*,*) '  -nea   write summary to CSV file in NEA benchmark format'
        write (*,*) '  -time  print timing summary'
        write (*,*) '  -tfu   print fuel temperature edits'
        write (*,*) '  -bwr   print BWR T/H edits'
        write (*,*) '  -batch print batch edits'
        write (*,*) '  -sM    edits a single statepoint M'
        write (*,*) '  -dN    distribution N (see below)'
        write (*,*)
        write (*,*) 'distribution selections for -dN option:'
        do idis=1, maxdist
          write (*,18) idis, trim(dist_label(idis))
        enddo
        stop
      endif
   18 format ('   -d',i0,5x,a)

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
        elseif (carg.eq.'-exit') then
          ifexit=.true.
        elseif (carg.eq.'-debug') then
          ifdebug=.true.
        elseif (carg.eq.'-csv') then
          ifcsv=.true.
        elseif (carg.eq.'-nea') then
          ifnea=.true.
        elseif (carg.eq.'-load') then   ! print loadings
          ifload=.true.
        elseif (carg.eq.'-time') then
          iftime=.true.
        elseif (carg.eq.'-tfu') then
          iftfu =.true.
        elseif (carg.eq.'-bwr') then
          ifbwr =.true.
        elseif (carg.eq.'-batch') then
          ifbatch=.true.
        elseif (carg.eq.'-help' .or. carg.eq.'--help') then  ! add for Ben
          write (*,*) 'run mpactread with no command line arguments for help'
        elseif (carg(1:2).eq.'-s') then   ! single statepoint output option
          read (carg(3:),*) istate
        elseif (carg(1:2).eq.'-d') then
          read (carg(3:),*,iostat=ierror) idis
          if (idis.ge.1 .and. idis.le.maxdist) dist_print(idis)=.true.
          if (ierror.ne.0) stop 'error reading distribution number on command line'
        else
          inputfile=carg
        endif
      enddo

      if (idis.eq.-1) then   ! print all distributions by default
        dist_print(:)=.true.
      endif

      if (iftfu) dist_print(lltfu)=.true.

      if (if1d)   write (*,*) '1D edits requested on command line'
      if (if2d)   write (*,*) '2D pin edits requested on command line'
      if (if2da)  write (*,*) '2DA assembly edits requested on command line'
      if (if2pin) write (*,*) '2PIN edits requested on command line'
      if (if3d)   write (*,*) '3D pin edits requested on command line'
      if (ifexit) write (*,*) 'exit maps requested on command line'
      if (istate.ne.0) write (*,*) 'single statepoint output selected'

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
        stop 'input file does not exist'
      endif
      call h5fopen_f (inputfile, H5F_ACC_RDONLY_F, file_id, ierror)   ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop 'H5 input file could not be opened'
      endif

!---------------------------
!  Read top level HDF data
!---------------------------

      xkeff=-100.0d0
      iver=-100

      idim(:)=0

! version

      dataset='veraout_version'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_integer(file_id, dataset, iver)
      else
        iver=0
        write (*,'(2a)') 'WARNING: veraout_version dataset does not exist on this file'
      endif

! title

      title=' '

      dataset='title'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call read_string(file_id, dataset, title)
      else
        write (*,'(2a)') 'WARNING: title dataset does not exist on this file'
      endif

      write (*,*)
      write (*,'(2a)')     ' Title: ', trim(title)
      write (*,'(a,i2)') ' VERAOUT version ', iver
      write (*,*)

!-------------------
!  Read CORE group
!-------------------

      call readcore(file_id, ifdebug)

!--- print loadings

      if (ifload) then
        dataset='Core Loading'
        if (if3d) then
          call print_3D_pin_map(dataset, npin, kd, nassm, pinload)
        endif
        if (if1d) then
!fix      call print1d(dataset, npin, kd, nassm, pinload, axial)   ! need to update routine to calculate loading
        endif

        if (if2d .or. if2da) then
          allocate (tdist2d(npin, npin, nassm))        ! 2D distributions

          iflag=.false.   ! do not weight with loading
          call collapse2d(dataset, npin, kd, nassm, axial, pinload, tdist2d, iflag)

          if (if2d) then
            call print_3D_pin_map('2D '//dataset, npin, 1,  nassm, tdist2d)
          endif
          if (if2da) then     ! 2D average assembly edits (already collapsed)
            call print2d_assm_map(dataset, npin, nassm, tdist2d)
          endif

          deallocate (tdist2d)
        endif

        if (if2pin) then
          call print2d_2pin_map(dataset, npin, nassm, pinload, kd)
        endif

      endif

!--- initialize batch edits

      if (ifbatch) then
        call batch_init(file_id, mapcore, icore, jcore)
      endif

!--- open NEA output files

      if (ifnea) then      ! NEA benchmark edits
        call nea_open()
      endif

!---------------------
!  Read STATE groups
!---------------------

!--- count number of statepoints

      if (istate.eq.0) then
        nstate=0
        do
          nstate=nstate+1
          write (group_name,'(a,i4.4,a)') '/STATE_', nstate, '/'
          if (ifdebug) write (*,*) 'debug: checking state= ', group_name

          call h5lexists_f(file_id, group_name, ifxst, ierror)
          if (.not.ifxst) then
            if (ifdebug) write (*,*) 'dataset not found - exiting'
            nstate=nstate-1
            exit
          endif
        enddo
      else
        nstate=1
      endif

      maxstate=nstate   ! actual number of statepoints to read

      allocate (state_xkeff(maxstate))
      allocate (state_expo(maxstate))  ! exposure
      allocate (state_efpd(maxstate))  ! exposure EFPD
      allocate (state_boron(maxstate))
      allocate (state_flow (maxstate))
      allocate (state_power(maxstate))
      allocate (state_tinlet(maxstate))
      allocate (state_3exp(maxstate))
      allocate (state_3pin(maxstate))
      allocate (state_2pin(maxstate))
      allocate (state_axoff(maxstate))
      allocate (state_nout(maxstate))
      allocate (state_time(maxstate))
      allocate (state_atfu(maxstate))
      allocate (state_3tfu(maxstate))

      state_xkeff(:)=0.0d0
      state_expo(:)=-100.0d0
      state_efpd(:)=-100.0d0
      state_boron(:)=-1.0d0
      state_flow (:)=-1.0d0
      state_power(:)=-1.0d0
      state_tinlet(:)=-1.0d0
      state_3exp(:)=0.0d0
      state_3pin(:)=0.0d0
      state_2pin(:)=0.0d0
      state_axoff(:)=0.0d0
      state_nout(:)=0
      state_time(:)=-1.0d0
      state_atfu(:)=0.0d0
      state_3tfu(:)=0.0d0


! 40 format (/,'--------------------------',&
!            /,'  Reading statepoint ', i0, &
!            /,'--------------------------')

      do nstate=1, maxstate

        if (istate.ne.0) then
          write (group_name,'(a,i4.4,a)') '/STATE_', istate, '/'
        else
          write (group_name,'(a,i4.4,a)') '/STATE_', nstate, '/'
        endif

        if (ifdebug) write (*,*) 'debug: state= ', group_name

!--- check if statepoint exists  (only necessary for istate case)

        call h5lexists_f(file_id, group_name, ifxst, ierror)
        if (.not.ifxst) then
          write (0,*)  'statepoint does not exist'
          exit
        endif

!!      write (*,40) nstate   ! statepoint is already printed in reading dataset edits

        write (*,*)

!---------------------------
!  Read statepoint scalars
!---------------------------

        xkeff=-1.0d0

        dataset=trim(group_name)//'exposure'
        call hdf5_read_double(file_id, dataset, state_expo(nstate))

        dataset=trim(group_name)//'exposure_efpd'
        call hdf5_read_double(file_id, dataset, state_efpd(nstate))

        dataset=trim(group_name)//'keff'
        call hdf5_read_double(file_id, dataset, xkeff)

        dataset=trim(group_name)//'boron'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, state_boron(nstate))
        endif

        dataset=trim(group_name)//'flow'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, state_flow(nstate))
          if (ifbwr) write (*,*) 'flow = ', state_flow(nstate)
        endif

        dataset=trim(group_name)//'power'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, state_power(nstate))
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

        dataset=trim(group_name)//'tinlet'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, state_tinlet(nstate))
          if (ifbwr) write (*,65) 'tinlet = ', state_tinlet(nstate),' C'
          if (ifbwr) write (*,65) 'tinlet = ', state_tinlet(nstate)*1.8d0+32.0d0,' F'
        endif
    65  format (1x,a,f12.4,a)

        dataset=trim(group_name)//'core_inlet_temp'       ! temp edit
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xtemp)
          if (ifbwr) write (*,65) 'core_inlet_temp  = ', xtemp,' C'
        endif

        dataset=trim(group_name)//'core_outlet_temp'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xtemp)
          if (ifbwr) write (*,65) 'core_outlet_temp = ', xtemp,' C'
        endif

        dataset=trim(group_name)//'core_inlet_density'       ! temp edit
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xtemp)
          if (ifbwr) write (*,*) 'core_inlet_density  = ', xtemp*0.001d0,' g/cc'
        endif

        dataset=trim(group_name)//'core_outlet_density'       ! temp edit
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xtemp)
          if (ifbwr) write (*,*) 'core_outlet_density = ', xtemp*0.001d0,' g/cc'
        endif

        dataset=trim(group_name)//'bypass'                    ! temp edit
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, xtemp)
          if (ifbwr) write (*,*) 'bypass = ', xtemp
        endif

        dataset=trim(group_name)//'pressure'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (ifxst) then
          call hdf5_read_double(file_id, dataset, pressure)
          if (ifbwr) then
            write (*,*) 'pressure = ', pressure,' MPa'
            write (*,*) 'pressure = ', pressure*145.03773800722d0,' psi'
          endif
        endif

        if (ifdebug) write (*,*) 'debug: keff = ', xkeff

!--- save statepoint values

        state_xkeff(nstate)=xkeff

        write (*,'(a, f10.4,a)') ' exposure =', state_expo(nstate),' GWD/MT'
        write (*,'(a, f10.4,a)') ' exposure =', state_efpd(nstate),' EFPD'
        write (*,'(a, f12.7)')   ' keff     =', xkeff

!--------------------------
!  Read statepoint arrays
!--------------------------

!--- allocate arrays if first statepoint

        if (nstate.eq.1) then
          if (ifdebug) then
            write (*,*) 'debug: allocating pin arrays'
            write (*,*) 'debug pin: nassm = ', nassm
            write (*,*) 'debug pin: kd    = ', kd
            write (*,*) 'debug pin: npin  = ', npin
          endif
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
                  tdist(i,j,k,n)=temp4d(n,k,i,j)    ! note i,j are not transposed
                enddo
              enddo
            enddo
          enddo

          if (idis.eq.llpow) then
            power=tdist        ! save 3D power - needed for mask
          endif

          if (idis.eq.lltfu) then    ! apply mask to fuel temperatures
            call masktfu(npin,  kd, nassm, power, tdist)
          endif

          call stat3d(dist_label(idis), npin,  kd, nassm, tdist, xave, xtemp, xtemp2)
          if (idis.eq.llpow) then
            state_3pin(nstate)=xtemp
            state_2pin(nstate)=xtemp2
            call calc_axoff (dist_label(idis), npin, kd, nassm, power, state_axoff(nstate))
          endif
          if (idis.eq.llexp) then
            state_3exp(nstate)=xtemp
          endif
          if (idis.eq.lltfu) then
            state_atfu(nstate)=xave   ! average
            state_3tfu(nstate)=xtemp  ! peak
          endif

          if (.not.dist_print(idis)) cycle    ! skip full edits if this was power

          if (if2d .or. if2da) then
            if (idis.eq.llexp) then
              iflag=.true.   ! collapse with pin loading
            else
              iflag=.false.
            endif
            call collapse2d(dist_label(idis), npin, kd, nassm, axial, tdist, tdist2d, iflag)
          endif

          if (idis.eq.llpow) then    ! special edits just for power
            call check_normalization(tdist)
          endif

          if (if3d) then
            call print_3D_pin_map(dist_label(idis), npin, kd, nassm, tdist)
          endif
          if (if2d) then
            call print_3D_pin_map('2D '//dist_label(idis), npin, 1,  nassm, tdist2d)
          endif
          if (if1d) then
            call print1d(dist_label(idis), npin, kd, nassm, tdist)
          endif

          if (if2da) then     ! 2D average assembly edits (already collapsed)
            call print2d_assm_map(dist_label(idis), npin, nassm, tdist2d)
          endif

!***** does 2pin use loading to collapse pin exposures?
          if (if2pin) then
            call print2d_2pin_map(dist_label(idis), npin, nassm, tdist, kd)
          endif

          if (ifexit) then
            call print_exit_map(dist_label(idis),  npin,  kd, nassm, tdist)
          endif

          if (ifbatch .and. idis.eq.llexp) then
             call batchstat (npin, kd, nassm, icore, jcore, mapcore, axial, tdist, power)
          endif

        enddo   ! distributions

        if (ifnea) then      ! NEA benchmark edits
          call print_nea(state_efpd(nstate), dist_label(llpow), nstate, npin,  kd, nassm, tdist)
        endif

        if (ifbatch) then    ! batch edits
          call batchedit
        endif

      enddo   ! end of statepoint loop - nstate

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------

      call h5fclose_f(file_id, ierror)

!--- close NEA output files

      if (ifnea) then      ! NEA benchmark edits
        call nea_close()
      endif

!--- deallocate distributions

      if (allocated(tdist2d)) then  ! protect from missing statepoints
        deallocate (tdist2d)
        deallocate (tdist)
        deallocate (power)
        deallocate (temp4d)
      endif

      deallocate (mapcore)
      deallocate (axial)

!--- print summary

      write (*,108)
      if (kd.eq.1) then   ! 2d
        write (*,112)
        do n=1, maxstate
          write (*,122) n, &
                 state_expo(n),  state_efpd(n), state_xkeff(n), &
                 state_boron(n), state_2pin(n), state_3exp(n)
        enddo
      else              ! 3D
        write (*,110)   ! 3D
        do n=1, maxstate
          write (*,120) n, &
                 state_expo(n),  state_efpd(n), state_xkeff(n), &
                 state_boron(n), state_2pin(n), state_3pin(n),  &
                 state_3exp(n),  state_axoff(n)
!x               state_flow(n), state_power(n), state_tinlet(n),
        enddo
      endif
  108 format (/,'==================================',&
              /,'       Statepoint Summary', &
              /,'==================================')
  110 format (  '   N   exposure  exposure  eigenvalue   boron      2PIN      3PIN      3EXP', &
                '     A/O(%)')
  120 format (i4, f10.4, f10.2, f12.6, f10.2, 7f10.4)

  112 format (  '   N   exposure  exposure  eigenvalue   boron      2PIN      2EXP')
  122 format (i4, f10.4, f10.2, f12.6, f10.2, 7f10.4)

  114 format (/,'   N   exposure  eigenvalue   power     flow     tinletC')
  124 format (i4, f10.4, f12.6,              7f10.4)

!--- BWR T/H edits

      if (ifbwr) then
        write (*,114)
        do n=1, maxstate
          write (*,124) n, &
                 state_expo(n),  state_xkeff(n), &
                 state_power(n), state_flow(n), state_tinlet(n)
        enddo
      endif

!--- fuel temperature edits

      if (iftfu) then
        write (*,220)
        do n=1, maxstate
          write (*,222) n, state_expo(n), state_efpd(n), state_atfu(n), state_3tfu(n)
        enddo
      endif
  220 format (/,'Fuel Temperature Edits (C)', &
              /,'   N   exposure  exposure    ave-tfu     max-tfu')
  222 format (i4, f10.4, f10.2, 2f12.4)

!--- print summary to CSV file

      if (ifcsv) then
        write (*,'(/,a)') ' writing output summary to "output.csv"'
        open (33,file='output.csv')
        write (33,'(a)') trim(inputfile)
        write (33,310)
        do n=1, maxstate
          write (33,320) n, &
                 state_expo(n),  state_efpd(n), state_xkeff(n), &
                 state_boron(n), state_2pin(n), state_3pin(n), state_3exp(n), &
                 state_axoff(n)
        enddo
        close (33)
      endif

  310 format ('Statepoint Summary', &
            /,'N, exposure, exposure, eigenvalue, boron, 2PIN, 3PIN, 3EXP, A/O(%)')
  320 format (i4,',', f10.4,',', f10.2,',', f12.6,',', f10.2,7(',',f10.4))

!--- timing summary

      if (iftime) then
        write (*,150)
        do n=1, maxstate
          write (*,160) n, state_expo(n), state_efpd(n), state_nout(n), state_time(n)
        enddo
      endif
  150 format (/,'==================================',&
              /,'       Timing Summary', &
              /,'==================================',&
              /,'   N   exposure  exposure  outers   time (sec)')
  160 format (i4, f10.4, f10.2, i8, 1x, f12.3)

!--- finished

      deallocate (state_xkeff)
      deallocate (state_expo)  ! exposure
      deallocate (state_efpd)  ! exposure EFPD
      deallocate (state_boron)
      deallocate (state_flow )
      deallocate (state_power)
      deallocate (state_tinlet)
      deallocate (state_3exp)
      deallocate (state_3pin)
      deallocate (state_2pin)
      deallocate (state_axoff)
      deallocate (state_nout)
      deallocate (state_time)
      deallocate (state_atfu)
      deallocate (state_3tfu)

      write (*,'(/,a)') 'done'

      end program

!=======================================================================
!
!  Debug subroutine to check power normalization
!
!=======================================================================
      subroutine check_normalization(power)
      use mod_coregeom, only : axial, icore, jcore, mapcore, npin, kd, nassm
      implicit none
      real(8), intent(in) :: power(npin,npin,kd,nassm)

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
!  MPACT has non-zero values for fuel temperatures in non-fuel rods
!  so zero the values out if the power is zero
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

