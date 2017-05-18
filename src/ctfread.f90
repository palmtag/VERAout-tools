   program ctfread
!=======================================================================
!
!  Program to read CTF HDF output file and print summary
!
!  Copyright (c) 2014-2017 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2014/05/16 - original version based on mpactread
!  2014/06/26 - updated code based on new CTF HDF format
!  2014/07/18 - updated power edits to use new W/cm^3 units
!               (this will likely change again)
!  2014/09/01 - change HDF open to read only
!  2014/08/25 - updated power edits to use new W/cm units (hopefully final time units change)
!  2014/09/22 - Update CTF data structure on HDF file
!  2014/10/10 - Added steaming rate and rod surface temperature arrays
!  2015/01/09 - Add multiple statepoint support
!  2015/03/06 - Major restructure to allow edits by distribution
!  2015/04/03 - Add DNB edits
!  2015/04/06 - Add core symmetry and core map
!  2015/12/09 - Update distributions
!  2017/04/01 - Update distributions
!
!  There are still some issues that need to be worked out:
!    * There are some cases that have a very small power in non-fuel regions
!      The user needs to examine the output carefully to make sure the
!      averages are correct.
!
!-----------------------------------------------------------------------
      use  hdf5
      use  mod_hdftools
      implicit none

      character(len=80)  :: filename
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i, j
      integer            :: idis
      integer            :: ierror
      integer            :: itype
      integer            :: ndim            ! number of dimensions in HDF file
      integer            :: idim(10)        ! value of dimensions on HDF file
      integer            :: nstate=0        ! statepoint number
      integer            :: lltcool
      integer            :: lltfu
      integer            :: llpow

      real(8)            :: prated
      real(8)            :: zsum
      real(8)            :: xtemp
      real(8)            :: xtemp2
      real(8)            :: qoi1, qoi2, qoi3, qoi4, qoi5

      logical            :: ifxst           ! flag if file or dataset exists
      logical            :: ifdebug=.false. ! debug flag

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=22)  :: state_name      ! HDF group name for STATE
      character(len=22)  :: group_name      ! HDF group name for CORE

      integer, parameter :: maxdist=8    ! was 19
      character(len=40)  :: dist_label(maxdist)
      logical            :: dist_print(maxdist)
      logical            :: dist_chan (maxdist)   ! channel array flag (i.e. not pin array)

      integer(hid_t)     :: file_id         ! HDF file ID number

! input data

      integer  :: isym                      ! core symmetry
      integer  :: icore, jcore              ! size of core map
      integer  :: kd                        ! number of axial levels in pin maps
      integer  :: nassm                     ! number of assemblies in pin maps
      integer  :: npin                      ! number of pins across one side of assembly
      integer  :: npin_save                 ! number of pins across one side of assembly
      integer  :: nchan                     ! number of channels across one side of assembly

      integer, allocatable :: mapcore(:,:)  ! core map (icore,jcore)
      real(8), allocatable :: axial(:)      ! axial elevations

      real(8), allocatable :: temp1(:,:,:,:)
      real(8), allocatable :: tdist(:,:,:,:)

      real(8), allocatable :: pow2d(:,:,:)
      real(8), allocatable :: charea(:,:,:,:)    ! channel areas

      real(8), allocatable :: power(:,:,:,:)
      real(8), allocatable :: tfuel(:,:,:,:)
      real(8), allocatable :: tcool(:,:,:,:)  ! coolant temperatures per pincell

! command line flags

      logical  :: if3d   =.false.           ! turn on 3D edits
      logical  :: if2d   =.false.           ! turn on 2D edits
      logical  :: if1d   =.false.           ! turn on 1D edits
      logical  :: iftfuel=.false.           ! turn on special fuel temp edits
      logical  :: ifexit =.false.           ! turn on channel exit edits

!  initialize

      i=0
      filename=' '
      dist_chan(:)=.false.   ! mark as pin arrays (not channel arrays)

      i=i+1
      dist_label(i) ="pin_powers"     ! no units
      llpow=i      ! save location of power
!old  i=i+1
!old  dist_label(i) ="Rod_Surface_Temp_NE_Quad [C]"   ! may only exist if steaming on
!old  i=i+1
!old  dist_label(i) ="Rod_Surface_Temp_NW_Quad [C]"   ! may only exist if steaming on
!old  i=i+1
!old  dist_label(i) ="Rod_Surface_Temp_SE_Quad [C]"
!old  i=i+1
!old  dist_label(i) ="Rod_Surface_Temp_SW_Quad [C]"
!old  i=i+1
!old  dist_label(i) ="Steaming_Rate_NE_Quad [kg_per_s]"
!old  i=i+1
!old  dist_label(i) ="Steaming_Rate_NW_Quad [kg_per_s]"
!old  i=i+1
!old  dist_label(i) ="Steaming_Rate_SE_Quad [kg_per_s]"
!old  i=i+1
!old  dist_label(i) ="Steaming_Rate_SW_Quad [kg_per_s]"
      i=i+1
      dist_label(i) ="pin_fueltemps [C]"
      lltfu=i      ! save location of tfuel
      i=i+1
      dist_label(i) ="pin_max_clad_surface_temp"
      i=i+1
      dist_label(i) ="pin_steamrate"
!old  i=i+1
!old  dist_label(i) ="pin_min_dnbr"
!old  i=i+1
!old  dist_label(i) ="equilibrium_quality"
!old  dist_chan(i) =.true.   ! mark as channel array
      i=i+1
      dist_label(i) ="channel_liquid_temps [C]"
      dist_chan(i) =.true.   ! mark as channel array
      lltcool=i    ! save location of tcool
      i=i+1
      dist_label(i) ="liquid_density"
      dist_chan(i) =.true.   ! mark as channel array
      i=i+1
      dist_label(i) ="mixture_massflux"
      dist_chan(i) =.true.   ! mark as channel array
      i=i+1
      dist_label(i) ="pressure"
      dist_chan(i) =.true.   ! mark as channel array
      i=i+1
      dist_label(i) ="vapor_void"
      dist_chan(i) =.true.   ! mark as channel array

      if (i.ne.maxdist) then
        write (*,*) 'number of distributions = ', i
        write (*,*) 'expecting ', maxdist
        stop 'maxdist error'
      endif

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  ctfread.exe [hdf5_file] {1D/2D/3D/tfuel/exit} {-dN}'

        write (*,*)
        write (*,*) 'list of distributions for -dN option:'
        do idis=1, maxdist
          write (*,18) idis, trim(dist_label(idis))
        enddo

        stop
      endif
   18 format ('   -d',i0,2x,a)

! parse command line arguments

      idis=-1
      dist_print(:)=.false.

      do i=1, iargs
        call get_command_argument(i,carg)
        if     (carg.eq.'3D' .or. carg.eq.'3d') then
          if3d=.true.
        elseif (carg.eq.'2D' .or. carg.eq.'2d') then
          if2d=.true.
        elseif (carg.eq.'1D' .or. carg.eq.'1d') then
          if1d=.true.
        elseif (carg.eq.'tfuel') then
          iftfuel=.true.
        elseif (carg.eq.'exit') then
          ifexit=.true.
        elseif (carg.eq.'debug') then
          ifdebug=.true.
        elseif (carg(1:2).eq.'-d') then
          read (carg(3:),*) idis
          if (idis.ge.1 .and. idis.le.maxdist) dist_print(idis)=.true.
        else
          filename=carg
        endif
      enddo

      if (idis.eq.-1) then   ! print all distributions if none chosen
        dist_print(:)=.true.
      endif

      if (iftfuel) then       ! turn on options needed for fuel temp fit
        dist_print(lltcool)=.true.   ! turn on coolant temp
        dist_print(lltfu)=.true.     ! turn on fuel temp
        dist_print(llpow)=.true.     ! turn on power
      endif


!--- initialize HDF fortran interface

      call h5open_f(ierror)

!--- open HDF file

      write (*,'(2a)') 'reading h5 file: ', trim(filename)
      inquire(file=filename, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(filename),' does not exist'
        stop
      endif
      call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, ierror)  ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(filename),' could not be opened'
        stop
      endif

!-------------------
!  Read CORE group
!-------------------

      group_name='/CORE/'
      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (.not.ifxst) then
        group_name=' '
        write (*,*) 'CORE group not found on HDF file - is this an old file?'
        stop 'CORE group not found on HDF file'
      endif

!--- symmetry

      dataset=trim(group_name)//'core_sym'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_integer(file_id, dataset, isym)
        write (*,*) 'core symmetry ', isym
      else
        isym=-1
      endif

!--- rated power

      dataset=trim(group_name)//'rated_power'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, prated)
        write (*,*) 'rated power ', prated
      else
        prated=0.0d0
      endif

!--- core map

      icore=0
      jcore=0

      dataset=trim(group_name)//'core_map'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
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
      endif

!--- read axial heights

      dataset=trim(group_name)//'channel_cell_height [cm]'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.1) stop 'invalid dimensions in cell heights'
      kd=idim(1)
      if (ifdebug) write (*,*) 'debug: kd = ', kd

      allocate (axial(kd))

      call hdf5_read_double(file_id, dataset, kd, i, axial)
      if (i.ne.kd) then
        stop 'invalid number of axial levels'
      endif

      if (ifdebug) then
        write (*,*) 'debug: axial heights [cm]'
        zsum=0.0d0
        do i=kd, 1, -1
           write (*,'(i4,f12.6)') i, axial(i)
           zsum=zsum+axial(i)
        enddo
        write (*,'(a,f10.5)') ' sum of axial heights =', zsum

      endif

!--- read channel flow areas (2D array)

      dataset=trim(group_name)//'channel_area'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.4) stop 'invalid dimensions in channel flow areas'

      write (*,*) 'channel area size ', idim(1:ndim)

      nassm=idim(1)    ! order (nassm, nchan, nchan)
      nchan=idim(3)
      npin=nchan-1

      if (ifdebug) then
        write (*,*) 'dimensions:'
        write (*,*) ' nassm = ', nassm
        write (*,*) ' npin  = ', npin
        write (*,*) ' nchan = ', nchan, idim(4)
        write (*,*) ' kd    = ', kd
      endif

      if (nassm.eq.0) stop 'invalid number of assemblies in pin data'
      if (npin .eq.0) stop 'invalid number of pins in pin data'
      if (idim(2).ne.kd)      stop 'invalid axial planes'
      if (idim(3).ne.idim(4)) stop 'invalid numbering npin by npin'

      allocate (charea(nchan,nchan,kd,nassm))    ! channel areas
      allocate (temp1(nassm,kd,nchan,nchan))
      call hdf5_read_double(file_id, dataset, nassm, kd, nchan, nchan, temp1)
      call transpose4d(nchan, nchan, kd, nassm, temp1, charea)
      deallocate (temp1)

!  check channel areas

      if (ifdebug) call check_channel(nchan, kd, nassm, charea, axial)

!--- allocate other arrays

      if (ifdebug) write (*,*) 'debug: allocating pin arrays'

      if (iftfuel) then
        allocate (power (npin, npin, kd, nassm))
        allocate (tfuel (npin, npin, kd, nassm))
        allocate (tcool (npin, npin, kd, nassm))  ! coolant temps per pincell
      endif

      npin_save=npin   ! save value for allocations below

!--- special check for old power label with units

      if (llpow.gt.0) then
        dataset='/STATE_0001/'//dist_label(llpow)
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,*) '*** Power label not found, checking old label ****'
          dist_label(llpow)="pin_powers [W per cm]"
          dataset='/STATE_0001/'//dist_label(llpow)
          call h5lexists_f(file_id, dataset, ifxst, ierror)
          if (.not.ifxst) stop 'no power label found'
        endif
      endif

!----------------------------------------------------
!  Read STATE group
!----------------------------------------------------

  40 format (/,'--------------------------',&
             /,'  Reading statepoint ', i0, &
             /,'--------------------------')

      nstate=0
      do
        nstate=nstate+1
        write (state_name,'(a,i4.4,a)') '/STATE_', nstate, '/'
        if (ifdebug) write (*,*) 'debug: state= ', state_name

!--- check if statepoint exists or data is in root group

        call h5lexists_f(file_id, state_name, ifxst, ierror)
        if (.not.ifxst) then
          if (ifdebug) write (*,*) 'dataset not found - exiting'
          if (nstate.eq.1) then   ! fatal error for first state only
            stop 'ERROR: First statepoint not found on file'
          endif
          goto 800
        endif

        write (*,40) nstate

!--------------------------------------------------------------------------------
! Read Statepoint Data
!--------------------------------------------------------------------------------

!--- read QOI (make separte input flag?)

      dataset=trim(state_name)//'QOI/cool_max_liquid_temp'
      call hdf5_read_double(file_id, dataset, qoi1)
      dataset=trim(state_name)//'QOI/core_avg_linear_heatrate'
      call hdf5_read_double(file_id, dataset, qoi2)
      dataset=trim(state_name)//'QOI/pin_max_clad_temp'
      call hdf5_read_double(file_id, dataset, qoi3)
      dataset=trim(state_name)//'QOI/pin_max_linear_power'
      call hdf5_read_double(file_id, dataset, qoi4)
      dataset=trim(state_name)//'QOI/pin_max_temp'
      call hdf5_read_double(file_id, dataset, qoi5)

      write (*,*)
      write (*,130) 'QOI cool_max_liquid_temp    ', qoi1, ' C'
      write (*,130) 'QOI core_avg_linear_heatrate', qoi2, ' W/cm'
      write (*,130) 'QOI pin_max_clad_temp       ', qoi3, ' C'
      write (*,130) 'QOI pin_max_linear_power    ', qoi4, ' W/cm'
      write (*,130) 'QOI pin_max_temp            ', qoi5, ' C'
      write (*,*)
  130 format (1x,a,f14.6,a)

!--- read distributions

      do idis=1, maxdist
        if (.not.dist_print(idis)) cycle   ! skip distribution

        dataset=trim(state_name)//dist_label(idis)
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') 'WARNING: dataset does not exist for this file ', trim(dataset)
          cycle
        endif

        if (dist_chan(idis)) then    ! set size of temporary array
           npin=nchan
        else
           npin=npin_save
        endif
        allocate (tdist(npin, npin, kd, nassm))

!--- read 4D data into temporary array then change order

        allocate (temp1(nassm,kd,npin,npin))
        dataset=trim(state_name)//dist_label(idis)
        call hdf5_read_double(file_id, dataset, nassm, kd, npin, npin, temp1)
        call transpose4d(npin, npin, kd, nassm, temp1, tdist)
        deallocate (temp1)

!--- check data

        if (idis.eq.llpow .or. idis.eq.lltfu) then
          call checkdata (dist_label(idis), npin, kd, nassm, tdist)
        endif

!--- save distributions for fuel temperature fit
!--- calculate pincell-based coolant temperatures from channel distributions

        if (iftfuel) then
          if (idis.eq.lltcool) then
            if (.not.dist_chan(idis)) stop 'failed sanity check'
            call pincell_coolant(nchan, npin_save, kd, nassm, charea, tdist, tcool)
          endif
          if (idis.eq.llpow) then
            power=tdist
          endif
          if (idis.eq.lltfu) then
            tfuel=tdist
          endif
        endif

!--- print overall statistics

        call stat3d_simple(dist_label(idis),  npin,  kd, nassm, axial, tdist, xtemp, xtemp2)
        write (*,*) '(coolant averages do not include flow area weighting)'

!--- 3D edits

        if (if3d) then
          call print_3D_pin_map(dist_label(idis),  npin,  kd, nassm, tdist)
        endif

!--- 2D edits

        if (if2d) then
          allocate (pow2d(npin, npin, nassm))
          call collapse2d(dist_label(idis), npin, kd, nassm, axial, tdist, pow2d)
          call print_3D_pin_map('2D '//dist_label(idis), npin, 1,  nassm, pow2d)
          deallocate (pow2d)
        endif

!--- print exit maps

        if (ifexit) then
          write (*,*)
          write (*,*) '===== Exit Maps ====='
          call print_exit_map(dist_label(idis),  npin,  kd, nassm, tdist)
        endif

!--- 1D edits

        if (if1d) then
          call print1d(dist_label(idis), npin, kd, nassm, tdist, axial)
        endif

        deallocate (tdist)

      enddo   ! loop over distributions

!=================================================================
!  Special edit to calculate fuel temperature fit
!=================================================================

      npin=npin_save     ! reset

!--- calculate quadratic fit of fuel temp vs. linear power

      if (iftfuel) then
        i=npin*npin*kd*nassm

!d      call quadratic(i, power, tfuel)  ! generate fit w/o subtracting coolant

!    calculate quadratic fit by subtracting coolant from fuel temp first
!    this edit also creates a large csv file

        call subtract_cool(npin, kd, nassm, power, tfuel, tcool)

        i=npin*npin*kd*nassm
        call quadratic(i, power, tfuel)

      endif

!--------------------------------------------------------------------------------
!   Finish Statepoints
!--------------------------------------------------------------------------------

      enddo     ! end of statepoint loop
  800 continue

!--- close hdf file

      call h5fclose_f(file_id, ierror)

!--- print summary here (?)

!   TO-DO:  save max values from each statepoint and print here

!--- deallocate memory

      deallocate (axial)
      deallocate (charea)

      if (iftfuel) then
        deallocate (tfuel)
        deallocate (tcool)
        deallocate (power)
      endif

      write (*,'(/,a)') 'done'

      end program

!=======================================================================
!
!  Transpose 4D array
!
!=======================================================================
      subroutine transpose4d(ni, nj, nk, nn, xold, xnew)
      implicit none
      integer, intent(in) :: ni, nj, nk, nn
      real(8) :: xold(nn, nk, nj, ni)
      real(8) :: xnew(ni, nj, nk, nn)

      integer :: i, j, k, n

      do n=1, nn
        do k=1, nk
          do j=1, nj
            do i=1, ni
              xnew(i,j,k,n)=xold(n,k,j,i)
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine transpose4d
!=======================================================================
      subroutine checkdata (label, npin, kd, nassm, power)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in) :: npin, kd, nassm
      real(8) :: power(npin,npin,kd,nassm)

      integer :: i, j, k, n
      integer :: nn   ! number of negative values
      integer :: nsm  ! number of small    values
      real(8) :: pp

      nn=0

      do n=1, nassm
        do k=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,n)
              if (pp.lt.0.0d0) then
                 nn=nn+1
!!               write (*,*) 'negative data ', label, i, j, k, n, pp
                 power(i,j,k,n)=0.0d0
              endif
              if (pp.gt.0.0d0 .and. pp.lt.1.0d-4) then
                 nsm=nsm+1
                 if (nsm.le.100) write (*,*) 'small data ', label, i, j, k, n, pp
                 power(i,j,k,n)=0.0d0
              endif
            enddo
          enddo
        enddo
      enddo

      if (nn.gt.0) then
        write (*,*) 'negative data found in ', nn, ' locations - fixing up'
        write (*,*) 'negative data is usually due to negative temperatures in guide tubes'
      endif
      if (nsm.gt.100) then
        write (*,*) 'small data found in ', nsm, ' locations - fixing up'
        write (*,*) 'only first 100 error messages printed'
      endif

      return
      end subroutine checkdata

!=======================================================================
!
!  Subroutine to calculate quadratic fit of fuel temperature vs. power
!
!    T(p) = c1 + c2 * p + c3 * p^2
!
!  See: http://www.personal.psu.edu/jhm/f90/lectures/lsq2.html
!
!=======================================================================
      subroutine quadratic(nsize, x, y)
      implicit none
      integer, intent(in) :: nsize
      real(8), intent(in) :: x(nsize)
      real(8), intent(in) :: y(nsize)

      integer :: nn, i

      real(8) :: sum1, sum2, sum3, sum4
      real(8) :: pp
      real(8) :: adet
      real(8) :: a(3,3)
      real(8) :: ainv(3,3)
      real(8) :: rhs(3), rnew(3)
      real(8) :: c(3)

      write (*,*)
      write (*,*) 'Generating quadratic fit for: T(p)=c1 + c2*p + c3*p*p'
      write (*,*) ' (does not include axial weighting)'

!--- calculate parameters only for non-zero power

      sum1=0.0d0
      sum2=0.0d0
      sum3=0.0d0
      sum4=0.0d0
      rhs(1)=0.0d0
      rhs(2)=0.0d0
      rhs(3)=0.0d0

      nn=0
      do i=1, nsize
        pp=x(i)
        if (pp.gt.0.0d0) then
           nn=nn+1
           sum1=sum1+pp
           sum2=sum2+pp*pp
           sum3=sum3+pp*pp*pp
           sum4=sum4+pp*pp*pp*pp
           rhs(1)=rhs(1)+y(i)
           rhs(2)=rhs(2)+y(i)*pp
           rhs(3)=rhs(3)+y(i)*pp*pp
        endif
      enddo
      if (nn.eq.0) then
        write (*,*) 'No power found in power arrays - is this a zero power case?'
        write (*,*) 'Skipping fuel temperature fit'
        return
      endif
      write (*,*) 'nn   = ', nn
      write (*,*) 'pave = ', sum1/dble(nn)
      write (*,*) 'tave = ', rhs(1)/dble(nn)

      sum1=sum1/dble(nn)
      sum2=sum2/dble(nn)
      sum3=sum3/dble(nn)
      sum4=sum4/dble(nn)
      rhs(:)=rhs(:)/dble(nn)

! setup:    A * c = rhs

      a(1,1) = 1.0d0
      a(1,2) = sum1
      a(1,3) = sum2
      a(2,1) = sum1
      a(2,2) = sum2
      a(2,3) = sum3
      a(3,1) = sum2
      a(3,2) = sum3
      a(3,3) = sum4

! invert 3x3 to get c()

      ainv = 0.0d+0

      ainv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      ainv(1,2) = A(3,2)*A(1,3) - A(1,2)*A(3,3)
      ainv(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)

      adet = ainv(1,1)*A(1,1) + ainv(1,2)*A(2,1) + ainv(1,3)*A(3,1)

      ainv(1,1) = ainv(1,1)/adet
      ainv(1,2) = ainv(1,2)/adet
      ainv(1,3) = ainv(1,3)/adet
      ainv(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/adet
      ainv(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/adet
      ainv(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/adet
      ainv(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/adet
      ainv(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/adet
      ainv(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/adet

! solve AINV*RHS to get coefficients

      c(1)=ainv(1,1)*rhs(1)+ainv(1,2)*rhs(2)+ainv(1,3)*rhs(3)
      c(2)=ainv(2,1)*rhs(1)+ainv(2,2)*rhs(2)+ainv(2,3)*rhs(3)
      c(3)=ainv(3,1)*rhs(1)+ainv(3,2)*rhs(2)+ainv(3,3)*rhs(3)

      write (*,*) 'c1 = ', c(1)
      write (*,*) 'c2 = ', c(2)
      write (*,*) 'c3 = ', c(3)

      write (*,*)
      write (*,*) 'Representative values of fit:'
      write (*,*) '  T(0)     =', c(1)
      write (*,*) '  T(pave)  =', c(1)+sum1*(c(2)+c(3)*sum1)
      write (*,*) '  T(pave*2)=', c(1)+sum1*2.0d0*(c(2)+c(3)*sum1*2.0d0)

! check matrix

      if (.false.) then

        rnew(1)=a(1,1)*c(1)+a(1,2)*c(2)+a(1,3)*c(3)
        rnew(2)=a(2,1)*c(1)+a(2,2)*c(2)+a(2,3)*c(3)
        rnew(3)=a(3,1)*c(1)+a(3,2)*c(2)+a(3,3)*c(3)

        write (*,*) 'check matrix:'
        write (*,*) '  old RHS and new RHS'
        write (*,*) rhs(1), rnew(1), rnew(1)-rhs(1)
        write (*,*) rhs(2), rnew(2), rnew(2)-rhs(2)
        write (*,*) rhs(3), rnew(3), rnew(3)-rhs(3)

      endif

      return
      end subroutine quadratic

!=======================================================================
      subroutine subtract_cool(npin, kd, nassm, power, tfuel, tcool)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8) :: power(npin, npin, kd, nassm)
      real(8) :: tfuel(npin, npin, kd, nassm)
      real(8) :: tcool(npin, npin, kd, nassm)

      integer :: i, j, na, k, m

      open (9,file='tfuel.csv')

      m=0
      do na=1, nassm
        do k=1, kd
          do j=1, npin
            do i=1, npin
              if (tfuel(i,j,k,na).gt.0.0d0) then
                tfuel(i,j,k,na)=tfuel(i,j,k,na)-tcool(i,j,k,na)
                m=m+1
                write (9,120) power(i,j,k,na), tfuel(i,j,k,na)
              endif
            enddo
          enddo
        enddo
      enddo


      write (*,*)
      write (*,*) 'Subtracting tcool from tfuel for ', m,' regions'

  120 format (f12.4,',', f12.4)

      close (9)
      return
      end subroutine subtract_cool

!=======================================================================
      subroutine check_channel(nchan, kd, nassm, charea, axial)
      implicit none
      integer, intent(in) :: nchan, kd, nassm
!old  real(8), intent(in) :: charea(nassm,nchan,nchan)
      real(8), intent(in) :: charea(nchan,nchan,kd,nassm)
      real(8), intent(in) :: axial(kd)

      integer :: i, j, na, k
      real(8) :: sum

      sum=0.0d0
      do k=1, kd
        sum=sum+axial(k)
      enddo

      write (*,*)
      write (*,*) 'Total core height = ', sum,' cm'

!--- edit flow areas

      do k=1, kd
      do na=1, nassm
        write (*,70) na
        do j=1, nchan
          write (*,'(20f7.4)') (charea(i,j,k,na), i=1, nchan)
        enddo

        sum=0.0d0
        do j=1, nchan
          do i=1, nchan
            sum=sum+charea(i,j,k,na)
          enddo
        enddo
        write (*,80) 'Total Channel flow area ', sum

      enddo
      enddo
   70 format (/,' Channel Flow Areas [cm^2]   - Assembly ',i0)
   80 format (1x,a,' =',f10.4,' [cm^2]')


!*** following is debug edit
!!    if (npin.eq.17) then
!!      apitch=21.54d0   ! needs to be from input
!!      write (*,80) 'Total assembly area         ', apitch**2 * nassm        ,'cm^2'
!!      write (*,80) 'flow should be approximately', apitch**2 * nassm * 0.5d0,'cm^2'
!!    endif

      write (*,*)

      return
      end subroutine check_channel

!=======================================================================
!
!  Create array of pincell-based coolant temperatures
!
      subroutine pincell_coolant(nchan, npin, kd, nassm, charea, chtemp, tcool)
      implicit none
      integer, intent(in) :: nchan, npin, kd, nassm
      real(8), intent(in) :: charea(nchan,nchan,kd,nassm)
      real(8), intent(in) :: chtemp (nchan,nchan,kd, nassm)  ! coolant temps per channel
      real(8), intent(out):: tcool(npin, npin, kd, nassm)    ! coolant temps per pincell

      real(8) :: parea(npin,npin,nassm)   ! automatic

      integer :: i, j, na, k
      real(8) :: sum
      real(8) :: w11, w12, w21, w22, ww, wsum

!--- calculate average channel temp per pincell  (needed for fuel temperature fit)

      sum=0.0d0
      wsum=0.0d0

      do j=1, npin
        do i=1, npin
          do na=1, nassm

            w11=1.0d0
            w12=1.0d0
            w21=1.0d0
            w22=1.0d0

            if (i.ne.1) then
              w11=w11*0.5d0
              w12=w12*0.5d0
            endif
            if (i.ne.npin) then
              w21=w21*0.5d0
              w22=w22*0.5d0
            endif

            if (j.ne.1) then
              w11=w11*0.5d0
              w21=w21*0.5d0
            endif
            if (j.ne.npin) then
              w12=w12*0.5d0
              w22=w22*0.5d0
            endif

            wsum=wsum+w11+w12+w21+w22

!!          write (*,'(2i3,a,2f8.4)') i, j, 'w=', w11, w12
!!          write (*,'(6x, a,4f8.4)')       '  ', w21, w22

         k=(kd+1)/2 !****** temp

            w11=w11*charea(i,j,    k,na)
            w12=w12*charea(i,j+1,  k,na)
            w21=w21*charea(i+1,j,  k,na)
            w22=w22*charea(i+1,j+1,k,na)

            ww=w11+w12+w21+w22
            sum=sum+ww         ! check sum
            parea(i,j,na)=ww   ! save for debug edit

            do k=1, kd

              tcool(i,j,k,na)=w11*chtemp(i,  j,  k,na) &
                              + w12*chtemp(i,  j+1,k,na) &
                              + w21*chtemp(i+1,j,  k,na) &
                              + w22*chtemp(i+1,j+1,k,na)
              tcool(i,j,k,na)=tcool(i,j,k,na)/ww

            enddo
          enddo
        enddo
      enddo

      write (*,*) 'Converting to pincell flow areas'
      write (*,80) 'Total Pincell flow area ', sum
   80 format (1x,a,' =',f10.4,' [cm^2]')

!d    do na=1, nassm
!d      write (*,*)
!d      write (*,*) 'Pincell Flow Areas [cm^2]   - Assembly ', na
!d      do j=1, npin
!d        write (*,'(20f7.4)') (parea(i,j,na), i=1, npin)
!d      enddo
!d    enddo

      return
      end subroutine pincell_coolant
!=======================================================================
!
!   Subroutine to convert 4 surface arrays into an array of max and average values
!
!   Do calculations inline so we don't have to declare any more arrays
!
!   input:
!     steam1   data for surface 1
!     steam2   data for surface 2
!     steam3   data for surface 3
!     steam4   data for surface 4
!
!   output:
!     steam1   max surface
!     steam2   average surface
!     steam3
!     steam4
!=======================================================================
      subroutine surfmaxave(nassm, npin, kd, steam1, steam2, steam3, steam4)
      implicit none
      integer, intent(in) :: nassm, npin, kd
      real(8)             :: steam1(npin, npin, kd, nassm)
      real(8)             :: steam2(npin, npin, kd, nassm)
      real(8)             :: steam3(npin, npin, kd, nassm)
      real(8)             :: steam4(npin, npin, kd, nassm)

      integer :: i, j, k, na
      real(8) :: temp1, temp2, temp3, temp4
      real(8) :: tmax     ! temporary max value
      real(8) :: tave     ! temporary ave value
      real(8) :: cmax     ! core wide max

      cmax=0.0d0
      do na=1, nassm
        do k=1, kd
          do j=1, npin
            do i=1, npin
              temp1=steam1(i,j,k,na)
              temp2=steam2(i,j,k,na)
              temp3=steam3(i,j,k,na)
              temp4=steam4(i,j,k,na)
              tmax=max(temp1,temp2,temp3,temp4)
              cmax=max(cmax,tmax)
              tave=(temp1+temp2+temp3+temp4)*0.25d0
              steam1(i,j,k,na)=tmax
              steam2(i,j,k,na)=tave
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine surfmaxave


