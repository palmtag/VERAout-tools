   module mod_coregeom
   implicit none
!=======================================================================
!
!  Module to read and store data from CORE block of HDF file
!
!  Copyright (c) 2014-2020 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================

      integer  :: isym                      ! core symmetry
      integer  :: icore, jcore              ! core size (assemblies)
      integer  :: nassm                     ! number of assemblies in core
      integer  :: kd
      integer  :: npin

      real(8)  :: rated_power               ! rated flow (kg/s)
      real(8)  :: rated_flow                ! rated power (MW)
      real(8)  :: apitch                    ! assembly pitch (cm)
      real(8)  :: coremass                  ! core_initial_mass (kg HM)
      real(8)  :: lhgr                      ! nominal_linear_heat_rate (W/cm)
      real(8)  :: axial0                    ! save bottom elevation (11.951 cm)

      integer, allocatable :: mapcore(:,:)  ! core assembly map
      integer, allocatable :: detector_map(:,:)  ! detector map (icore,jcore)

      real(8), allocatable :: axial(:)      ! axial elevations
      real(8), allocatable :: pinload(:,:,:,:)  ! pin loadings (npin,npin,kd,nassm)

      character(len=2), allocatable :: xlabel(:)  ! assembly map labels
      character(len=2), allocatable :: ylabel(:)  ! assembly map labels

   contains

!=======================================================================
!
!   Subroutine to read CORE block from VERAout HDF file
!
!=======================================================================
      subroutine readcore(file_id, ifdebug)
      use hdf5
      use mod_hdftools
      implicit none

      integer(hid_t), intent(in) :: file_id
      logical       , intent(in) :: ifdebug

!--- local

      integer :: i, j, k, n
      integer :: ierror
      integer :: itype
      integer :: naxial          ! number of axial levels
      integer :: ndim            ! temp variable
      integer :: idim(10)        ! temp variable

      logical :: ifxst

      real(8) :: rf2             ! flow in english units

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=12)  :: group_name

      real(8), allocatable :: temp4d(:,:,:,:)

!--- initialize

      isym=-100           ! core symmetry
      naxial=0            ! number of axial edit bounds
      icore=0             ! number of assemblies across
      jcore=0             ! number of assemblies across

      apitch=21.5d0       ! set default value
      rated_power=0.0d0   ! set default value
      rated_flow =0.0d0   ! set default value
      coremass   =0.0d0   ! set default value
      lhgr       =0.0d0   ! set default value

      idim(:)=0

!-------------------
!  Read CORE group
!-------------------

! check if CORE group exists (group may not exist on old MPACT files and all data will be in root)

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

!--- mass

      dataset=trim(group_name)//'core_initial_mass'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, coremass)
      else
        write (*,*) '**** core_initial_mass is missing from file ****'
      endif

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

      dataset=trim(group_name)//'nominal_linear_heat_rate'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, lhgr)
      else
        write (*,*) '**** nominal_linear_heat_rate is missing from file ****'
      endif

      rf2=rated_flow*7936.641438d-6   ! convert to lbm/hr

      write (*,*)
      if (isym.eq.1) then
        write (*,'(a)') ' core symmetry  isym 1 (full)'
      elseif (isym.eq.4) then
        write (*,'(a)') ' core symmetry  isym 4 (qtr)'
      else
        write (*,'(a,i2,a)') ' core symmetry  isym  = ', isym,' (unknown)'
      endif
      write (*,30)       'core mass      ', coremass,   ' kg HM'
      write (*,30)       'rated power    ', rated_power,' MW'
      write (*,30)       'rated flow     ', rated_flow, ' kg/s'
      write (*,30)       '               ', rf2,        ' Mlb/hr'
      write (*,30)       'nominal LHGR   ', lhgr,       ' W/cm'
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

!--- detector map

      allocate (detector_map(icore,jcore))
      detector_map=0

      dataset=trim(group_name)//'detector_map'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.2) stop 'invalid dimensions in detector_map'
      call hdf5_read_integer(file_id, dataset, icore, jcore, detector_map)

!d    write (*,'(/,a)') ' Detector Map:'
!d    do j=1, jcore
!d       write (*,'(2x,20i3)') (detector_map(i,j),i=1,icore)
!d    enddo

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
      do j=naxial+1, 1, -1
        write (*,'(i5,20f12.4)') j-1, axial(j)
      enddo
      if (naxial.eq.0) then
        write (*,*) '*** warning: no axial edits found ****'
      endif

      write (*,'(a,f12.4)') ' total core height = ', axial(naxial+1)-axial(1)
      write (*,*)

      axial0=axial(1)    ! save bottom height
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

!--- pin loadings

      dataset=trim(group_name)//'initial_mass'
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

      allocate (temp4d(nassm, kd, npin, npin))   ! temporary
      allocate (pinload(npin, npin, kd, nassm))

      call hdf5_read_double(file_id, dataset, nassm, kd, npin, npin, temp4d)

      do n=1, nassm
        do k=1, kd
          do j=1, npin
            do i=1, npin
              pinload(i,j,k,n)=temp4d(n,k,i,j)  ! note i,j are not transposed
            enddo
          enddo
        enddo
      enddo

      deallocate (temp4d)

      call fixload()       ! fix pin loading

      call coremass_pinsum()  ! find total core mass from pin loadings

!--- temp define edit labels

  !  xlabel  R P N M L K J H G  F  E  D  C  B  A
  !  ylabel  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

      allocate (xlabel(icore))
      allocate (ylabel(jcore))

!**** still having issues with reading strings - use default values instead

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

      end subroutine readcore

!=======================================================================
!
!  Fix error in MPACT pin loadings
!
!  The loadings are not expanded to full-assembly on the edge assemblies
!  like the other pin exposures.  This causes a lot of problems when
!  trying to do exposure collapses
!
!=======================================================================

      subroutine fixload()
      implicit none

      integer :: i, j, k
      integer :: ia, ja, na
      integer :: iac, jac
!     integer :: nl

      if (isym.ne.4) return   ! only do for qtr-core

      if (mod(icore,2).eq.0) return   ! skip for even cores

      iac=(icore+1)/2   ! center
      jac=(jcore+1)/2

!--- fix top edge

      ja=jac  ! top row, fold up
      do ia=iac, icore
        na=mapcore(ia,ja)
        do k=1, kd
          do j=1, (npin+1)/2
            do i=1, npin
              if (pinload(i,j,k,na).eq.0.0d0) then
                pinload(i,j,k,na)=pinload(i,npin-j+1,k,na)
              endif
            enddo
          enddo
        enddo
      enddo

!--- fix left edge

      ia=iac  ! top row, fold up
      do ja=jac, jcore
        na=mapcore(ia,ja)
        do k=1, kd
          do j=1, npin
            do i=1, (npin+1)/2
               if (pinload(i,j,k,na).eq.0.0d0) then
                 pinload(i,j,k,na)=pinload(npin-i+1,j,k,na)
               endif
            enddo
          enddo
        enddo
      enddo

!--- double check final results

!     write (*,*) 'final check on loading:'
!     do na=1, nassm
!       do k=1, kd
!         nl=0
!         do j=1, npin
!           do i=1, npin
!              if (pinload(i,j,k,na).gt.0.0d0) nl=nl+1
!           enddo
!         enddo
!         if (nl.ne.264) then
!           write (*,*) 'ERROR: final loading error na, k, nload ', na, k, nl
!         endif
!       enddo
!     enddo

      return
      end subroutine fixload
!=======================================================================
!
!  Subroutine to calculate total core mass from pin loading
!   ** no longer needed because total mass is now written to HDF file
!
!=======================================================================
      subroutine coremass_pinsum
      implicit none

      integer :: i, j, k
      integer :: ia, ja, na
      real(8) :: totk    ! total in plane
      real(8) :: total

      total=0.0d0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
          na=mapcore(ia,ja)
          if (na.eq.0) cycle

          do k=1, kd
            totk=0.0d0
            do j=1, npin
              do i=1, npin
                totk=totk+pinload(i,j,k, na)
              enddo
            enddo
            total=total+totk*axial(k)
          enddo

        enddo  ! ia
      enddo    ! ja

      if (total.gt.1000.0d0) then
        write (*,110) total
      else
        write (*,112) total
      endif
  110 format (' Total core mass ', f12.3,' kg HM')
  112 format (' Total core mass ', 1p,e14.5,' kg HM')

      return
      end subroutine coremass_pinsum
!=======================================================================
   end module mod_coregeom
