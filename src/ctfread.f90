!=======================================================================
!
!  Program to read CTF HDF output file and print summary
!
!  Copyright (c) 2014 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2014/05/16 - original version based on mpactread
!
!  2014/06/26 - updated code based on new CTF HDF format
!
!  2014/07/18 - updated power edits to use new W/cm^3 units
!               (this will likely change again)
!
!  2014/09/01 - change HDF open to read only
!
!  2014/08/25 - updated power edits to use new W/cm units (hopefully final time units change)
!
!  2014/09/19 - CTF output file changed again!  now the top level is "Simulation Results"
!
!  There are still some issues that need to be worked out:
!    * There are some cases that have a very small power in non-fuel regions
!      The user needs to examine the output carefully to make sure the 
!      averages are correct.
!
!-----------------------------------------------------------------------
      program ctfread
      use  hdf5
      use  mod_hdftools
      implicit none

      character(len=80)  :: filename
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i
      integer            :: ierror
      integer            :: itype
      integer            :: ndim            ! number of dimensions in HDF file
      integer            :: idim(10)        ! value of dimensions on HDF file
      integer            :: nstate=0        ! statepoint number

      real(8)            :: zsum
      real(8), parameter :: t_btu_J = 1055.05585262d0
      real(8), parameter :: cvert =t_btu_J/(3600.0d0*12.0d0*2.54d0)  ! W/cm to BTU/hr/ft

      logical            :: ifxst
      logical            :: ifdebug=.true.  ! debug flag

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=22)  :: state_name      ! HDF group name for STATE
      character(len=22)  :: group_name      ! HDF group name for CORE
      character(len=30)  :: label_power     ! power label - for both input and output
      character(len=30)  :: label_tfuel     ! tfuel label - for both input and output

      integer(hid_t)     :: file_id         ! HDF file ID number

! input data

      integer  :: kd                        ! number of axial levels in pin maps
      integer  :: nassm                     ! number of assemblies in pin maps
      integer  :: npin                      ! number of pins across one side of assembly
      integer  :: nchan                     ! number of channels across one side of assembly

      real(8), allocatable :: axial(:)      ! axial elevations

      real(8), allocatable :: temppower(:,:,:,:)

      real(8), allocatable :: power(:,:,:,:)
      real(8), allocatable :: tfuel(:,:,:,:)
      real(8), allocatable :: charea(:,:,:)    ! channel areas
      real(8), allocatable :: chtemp(:,:,:,:)  ! channel heights
      real(8), allocatable :: tcool (:,:,:,:)  ! coolant temperatures per pincell

! command line flags

      logical  :: if3d   =.false.           ! turn on 3D edits
      logical  :: if2d   =.false.           ! turn on 2D edits
      logical  :: if1d   =.false.           ! turn on 1D edits
      logical  :: iftfuel=.false.           ! turn on special fuel temp edits

!  initialize

      filename=' '

      label_power='pin_powers [W per cm]'
      label_tfuel='pin_fueltemps [C]'

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  ctfread.exe [hdf5_file] {1D/2D/3D/tfuel}'
        stop
      endif

! parse command line arguments

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
        else
          filename=carg
        endif
      enddo

      if (.not.if1d)  write (*,*) 'no 1D edits requested on command line'
      if (.not.if2d)  write (*,*) 'no 2D edits requested on command line'
      if (.not.if3d)  write (*,*) 'no 3D edits requested on command line'

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
        write (*,'(3a)') 'error: H5 input file ',trim(filename), &
                   ' could not be opened'
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

!--- read channel flow areas (3D array)

      dataset=trim(group_name)//'channel_flow_areas [cm^2]'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.3) stop 'invalid dimensions in channel flow areas'

      nassm=idim(1)    ! order (nassm, kd, npin, npin)
      nchan=idim(2)
      npin=nchan-1

      if (ifdebug) then
        write (*,*) 'dimensions:'
        write (*,*) ' nassm = ', nassm
        write (*,*) ' npin  = ', npin
        write (*,*) ' nchan = ', nchan, idim(3)
        write (*,*) ' kd    = ', kd   
      endif

      if (nassm.eq.0) stop 'invalid number of assemblies in pin data'
      if (npin .eq.0) stop 'invalid number of pins in pin data'
      if (idim(2).ne.idim(3)) stop 'invalid numbering npin by npin'

      allocate (charea(nassm,nchan,nchan))    ! channel areas
      call hdf5_read_double(file_id, dataset, nassm, nchan, nchan, charea)

!--- allocate other arrays

      if (ifdebug) write (*,*) 'debug: allocating pin arrays'

      allocate (power    (npin, npin, kd, nassm))
      allocate (tfuel    (npin, npin, kd, nassm))
      allocate (chtemp   (nchan,nchan,kd, nassm))
      allocate (tcool    (npin, npin, kd, nassm))  ! coolant temps per pincell

!----------------------------------------------------
!  Read STATE group - only one statepoint supported
!----------------------------------------------------

      nstate=1
      write (state_name,'(a,i4.4,a)') '/STATE_', nstate, '/'
      write (state_name,'(a,i4.4,a)') '/Simulation Results/'   ! ********

      if (ifdebug) write (*,*) 'debug: state= ', state_name

!--- check if statepoint exists or data is in root group

      call h5lexists_f(file_id, state_name, ifxst, ierror)
      if (.not.ifxst) then
!x      write (*,*) 'dataset not found - exiting'
!x      goto 800
        write (*,'(3a)') 'WARNING: Group ', trim(state_name),' not found'
        write (*,'(3a)') 'WARNING: Looking in HDF root instead'
        state_name=' '
      endif

!--------------------------------------------------------------------------------
! Read pin powers
!--------------------------------------------------------------------------------

!--- read 4D data into temporary array then change order

        allocate (temppower(nassm,kd,npin,npin))

        dataset=trim(state_name)//label_power
        call hdf5_read_double(file_id, dataset, nassm, kd, npin, npin, temppower)
        call transpose4d(npin, npin, kd, nassm, temppower, power)
        call checkdata ('power', npin, kd, nassm, power)

        dataset=trim(state_name)//label_tfuel
        call hdf5_read_double(file_id, dataset, nassm, kd, npin, npin, temppower)
        call transpose4d(npin, npin, kd, nassm, temppower, tfuel)
        call checkdata ('tfuel', npin, kd, nassm, tfuel)

        deallocate (temppower)

        allocate (temppower(nassm,kd,nchan,nchan))

        dataset=trim(state_name)//'channel_liquid_temps [C]'
        call hdf5_read_double(file_id, dataset, nassm, kd, nchan, nchan, temppower)
        call transpose4d(nchan, nchan, kd, nassm, temppower, chtemp)

        deallocate (temppower)

!  check channel areas

        call check_channel(nchan, npin, kd, nassm, charea, axial, chtemp, tcool)

!------------------
!    Edits
!------------------

!--- print 3D maps

        if (if3d) then
          call print_pin_map(label_power,                     npin,  kd, power,  nassm)
          call print_pin_map(label_tfuel,                     npin,  kd, tfuel,  nassm)
          call print_pin_map('Channel Flow Temperatures [C]', nchan, kd, chtemp, nassm)
        endif

        call stat3d(label_power, npin, kd, nassm, axial, power)        
        call stat3d(label_tfuel, npin, kd, nassm, axial, tfuel)

        call stat3d('Pincell Coolant Temperatures [C]', npin, kd, nassm, axial, tcool)
        write (*,*) '(coolant averages do not include flow area weighting)'

        call stat3d('Channel Coolant Temperatures [C]', nchan, kd, nassm, axial, chtemp)
        write (*,*) '(coolant averages do not include flow area weighting)'

!--- 2D edits

   ! **** To-Do

!--- 1D edits

        if (if1d) then
           call print1d(label_power, npin, kd, nassm, power, axial)
           call print1d(label_tfuel, npin, kd, nassm, tfuel, axial)
        endif

!=================================================================
!  Special edit to calculate fuel temperature fit
!=================================================================

!--- calculate quadratic fit of fuel temp vs. linear power

        if (iftfuel) then
          i=npin*npin*kd*nassm

!!        call quadratic(i, power, tfuel)  ! generate fit w/o subtracting coolant

!    calculate quadratic fit by subtracting coolant from fuel temp first
!    this edit also creates a large csv file

          call subtract_cool(npin, kd, nassm, power, tfuel, tcool)

          i=npin*npin*kd*nassm
          call quadratic(i, power, tfuel)

        endif

!--------------------------------------------------------------------------------
!   Finish Statepoints
!--------------------------------------------------------------------------------
!x800 continue

!--- close hdf file

      call h5fclose_f(file_id, ierror)

!--- deallocate memory

      if (allocated(tfuel)) then
        deallocate (tfuel)
        deallocate (power)
        deallocate (axial)
        deallocate (charea)
        deallocate (chtemp)
      endif

      write (*,'(/,a)') 'done'

      end program

!=======================================================================
!
!  Print 3D statistics
!
!=======================================================================
      subroutine stat3d (title, npin, kd, nassm, axial, power)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer            :: i, j, k
      integer            :: ia, kpin
      real(8)            :: pp
      real(8)            :: zave, zlen, zmin, zmax   ! 3D values

!--- calculate 3D statistics to set the format correctly

      kpin=0
      zave=0.0d0
      zlen=0.0d0
      zmin=1.0d20
      zmax=0.0d0

      do k=1, kd
        do ia=1, nassm
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,ia)
              if (pp.gt.0.0d0) then
                zave=zave+axial(k)*pp
                zlen=zlen+axial(k)
                kpin=kpin+1
                zmin=min(zmin,pp)
                zmax=max(zmax,pp)
              endif
            enddo
          enddo
        enddo
      enddo
      if (zlen.gt.0.0d0) then
        zave=zave/zlen
      else
        zmin=0.0d0    ! avoid overflow
      endif

      write (*,'(/,1x,2a)') '3D Statistics - ', trim(title)
      write (*,'(a,2f12.4)') ' Max ', zmax
      write (*,'(a,2f12.4)') ' Ave ', zave
      write (*,'(a,2f12.4)') ' Min ', zmin
      write (*,'(a,i12)  ')  ' Num ', kpin
!!    write (*,*) 'debug: zmin = ', zmin   

!--- check for very small values

      do k=1, kd
        do ia=1, nassm
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,ia)
              if (pp.gt.0.0d0) then
                if (pp.lt.zave*1.0d-4) then
                   write (*,*) 'small value detected ', i, j, k, ia, pp
                endif
              endif
            enddo
          enddo
        enddo
      enddo


      return
      end subroutine stat3d

!=======================================================================
!
!  Print pin powers by assembly using pretty output
!
!=======================================================================
      subroutine print_pin_map(title, npin, kd, power, nassm)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer            :: i, j, m
      integer            :: ia, klev, kpin
      real(8)            :: pp, pave, pmin, pmax
      real(8)            :: zmax   ! 3D max
      character(len=150) :: line
      character(len=8)   :: fmt

!--- calculate 3D statistics to set the format correctly
!---  **** average does not use axial weighting

      zmax=0.0d0

      do ia=1, nassm
        do klev=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,klev,ia)
              zmax=max(zmax,pp)
            enddo
          enddo
        enddo
      enddo

      write (*,'(/,1x,a)') trim(title)
!!    write (*,'(a,2f12.4)') ' Max ', zmax

!--- write maps

      fmt='(f7.4)'
      if (zmax.gt.100.0d0) fmt='(f7.1)'

      if (npin.gt.20) then
        write (*,*) '**** too many pins across to print nice maps ****'
        return
      endif

!**** To-do: read labels from input and use in edits


      do ia=1, nassm        ! loop over assemblies
        do klev=kd, 1, -1   ! loop over axial levels
          write (*,'(/,1x,a)') trim(title)
!!        write (*,'(1x,a,i3,2a)') 'Assembly ', ia,'    type ', trim(assmmap(nassm))   ! no types to mpact
          write (*,'(1x,a,i3,2a)') 'Assembly ', ia
          if (kd.gt.1) then
             write (*,'(  1x,a,i3)') 'Level    ', klev
          else
             write (*,*) '2D collapse'
          endif
          pave=0.0
          kpin=0
          pmin=2.0d20
          pmax=0.0d0
          do i=1, npin
            line=' '
            m=0
            do j=1, npin
              pp=power(i,j,klev,ia)
              if (pp.eq.0.0) then
                line(m+1:m+7)='   --- '
              else
                write (line(m+1:m+7),fmt) pp
                pave=pave+pp
                kpin=kpin+1
                pmin=min(pmin,pp)
                pmax=max(pmax,pp)
              endif
              m=m+7
            enddo
            write (*,'(1x,a)') line(1:m)
          enddo
          if (kpin.gt.0) pave=pave/real(kpin)
          if (pmin.gt.1.0d20) pmin=0.0d0   ! protect if no data
          write (*,210) kpin, pave, pmin, pmax
          write (*,*) 'debug average = ', pave

        enddo    ! klev
      enddo      ! ia
! 210 format (5x,'number of hot pins',i4,'    average=',f7.4,'   min=', f7.4,'   max=', f7.4)
  210 format (5x,'number of hot pins',i4,'    average=',f12.4,'   min=', f12.4,'   max=', f12.4)

      return
      end subroutine print_pin_map
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
                 write (*,*) 'small data ', label, i, j, k, n, pp
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
      write (*,*) 'T(0)     =', c(1)
      write (*,*) 'T(pave)  =', c(1)+sum1*(c(2)+c(3)*sum1)
      write (*,*) 'T(pave*2)=', c(1)+sum1*2.0d0*(c(2)+c(3)*sum1*2.0d0)

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
      subroutine check_channel(nchan, npin, kd, nassm, charea, axial, chtemp, tcool)
      implicit none
      integer, intent(in) :: nchan, npin, kd, nassm
      real(8), intent(in) :: charea(nassm,nchan,nchan)
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: chtemp (nchan,nchan,kd, nassm)  ! coolant temps per channel
      real(8), intent(out):: tcool(npin, npin, kd, nassm)  ! coolant temps per pincell

      real(8) :: parea(npin,npin,nassm)   ! automatic

      integer :: i, j, na, k
      real(8) :: sum
      real(8) :: w11, w12, w21, w22, ww, wsum

      sum=0.0d0
      do k=1, kd
        sum=sum+axial(k)
      enddo

      write (*,*)
      write (*,*) 'Total core height = ', sum,' cm'

!--- edit flow areas

      do na=1, nassm
        write (*,70) na
        do j=1, nchan
          write (*,'(20f7.4)') (charea(na,i,j), i=1, nchan)
        enddo

        sum=0.0d0
        do j=1, nchan
          do i=1, nchan
            sum=sum+charea(na,i,j)
          enddo
        enddo 
        write (*,80) 'Total Channel flow area ', sum

      enddo
   70 format (/,' Channel Flow Areas [cm^2]   - Assembly ',i0)
   80 format (1x,a,' =',f10.4,' [cm^2]')


!*** following is debug edit
!!    if (npin.eq.17) then
!!      apitch=21.54d0   ! needs to be from input
!!      write (*,80) 'Total assembly area         ', apitch**2 * nassm        ,'cm^2'
!!      write (*,80) 'flow should be approximately', apitch**2 * nassm * 0.5d0,'cm^2'
!!    endif


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

            w11=w11*charea(na,i,j)
            w12=w12*charea(na,i,j+1)
            w21=w21*charea(na,i+1,j)
            w22=w22*charea(na,i+1,j+1)

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

!d    do na=1, nassm
!d      write (*,*)
!d      write (*,*) 'Pincell Flow Areas [cm^2]   - Assembly ', na
!d      do j=1, npin 
!d        write (*,'(20f7.4)') (parea(i,j,na), i=1, npin)
!d      enddo
!d    enddo

      return
      end subroutine check_channel
!=======================================================================
!
!  Print axial edits
!
!=======================================================================
      subroutine print1d(title, npin, kd, nassm, power, axial)
      implicit none
      character(len=*), intent(in) :: title
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)

      integer          :: i, j, k, nn
      integer          :: na
      real(8)          :: z1, z2
      real(8)          :: pp
      real(8)          :: pave
      real(8)          :: zave
      real(8)          :: ztot, zheat
      real(8)          :: axpow(kd)    ! automatic
      real(8)          :: axmid(kd)    ! automatic
      real(8)          :: axmin(kd)    ! automatic
      real(8)          :: axmax(kd)    ! automatic
      integer          :: numax(kd)    ! automatic

      if (kd.le.1) return   ! skip 1D edits

!--- calculate axial averages **** this will be wrong if qtr-core and assemblies repeated

      axpow(:)=0.0d0
      axmin(:)=1.0d20
      axmax(:)=0.0d0
      numax(:)=0

      do na=1, nassm
        do k=1, kd       ! loop over axial levels
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                axpow(k)=axpow(k)+pp
                numax(k)=numax(k)+1
                axmin(k)=min(pp,axmin(k))  ! don't count zero power locations
                axmax(k)=max(pp,axmax(k))
              endif
            enddo
          enddo
        enddo    ! k
      enddo

!  Normalize each axial level
!  Note: this will give "odd" results if different axial levels have different number of pins

      pave=0.0d0   ! core average
      zave=0.0d0   ! height of fuel
      nn=0
      do k=1, kd
        if (numax(k).eq.0) then
           axmin(k)=0    ! protect overflow
        else
           axpow(k)=axpow(k)/dble(numax(k))
           pave=pave+axial(k)*axpow(k)
           zave=zave+axial(k)
           if (nn.eq.0) nn=numax(k)       ! save first non-zero level
           if (nn.ne.numax(k)) nn=-1000   ! set error
        endif
      enddo
      if (zave.gt.0.0d0) pave=pave/zave
      if (nn.lt.0) then
        write (*,*) '*** 1D axial distributions not calculated because some ***'
        write (*,*) '***    axial levels have different numbers of pins     ***'
        return
      endif

!--- calculate heated length

      ztot=0.0d0
      zheat=0.0d0
      do k=1, kd
        ztot=ztot+axial(k)
        if (axpow(k).gt.0.0d0) zheat=zheat+axial(k)
      enddo

!--- find axial midpoints for plotting

      z1=0.0d0
      do k=1, kd
        z2=z1+axial(k)
        axmid(k)=0.5d0*(z1+z2)
        z1=z2
      enddo

      write (*,'(/,1x,a,/,1x,a)') trim(title), '1D Axial Distribution'
      write (*,'(a)')      '    K   Elev         Delta      Average     Minimum     Maximum'
      do k=kd, 1, -1
        write (*,200) k, axmid(k), axial(k), axpow(k), axmin(k), axmax(k)
      enddo
      write (*,205) pave

      write (*,140) 'total ', ztot
      write (*,140) 'heated', zheat
  140 format (2x,a,' axial length ', f12.4,' cm')

  200 format (i5,f10.4,6f12.4)
  205 format (2x,'ave',22x,f12.4)

      return
      end subroutine print1d
