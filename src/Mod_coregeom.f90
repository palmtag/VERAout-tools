   module mod_coregeom
   implicit none
!=======================================================================
!
!  Module to read and store data from CORE block of HDF file
!
!  Copyright (c) 2014-2017 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================

      integer  :: isym                      ! core symmetry
      integer  :: naxial                    ! number of axial levels
      integer  :: icore, jcore              ! core size

      real(8)  :: rated_power               ! rated flow
      real(8)  :: rated_flow                ! rated power
      real(8)  :: apitch                    ! assembly pitch

      integer, allocatable :: mapcore(:,:)  ! core map

      real(8), allocatable :: axial(:)      ! axial elevations

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

      integer :: i, j, n
      integer :: ierror
      integer :: itype
      integer :: ndim            ! temp variable
      integer :: idim(10)        ! temp variable

      logical :: ifxst

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=12)  :: group_name

!-------------------
!  Read CORE group
!-------------------

      isym=-100           ! core symmetry
      naxial=0            ! number of axial edit bounds
      icore=0             ! number of assemblies across
      jcore=0             ! number of assemblies across

      apitch=21.5d0       ! set default value
      rated_power=0.0d0   ! set default value
      rated_flow =0.0d0   ! set default value

      idim(:)=0

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
      do j=naxial+1, 1, -1
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

      end subroutine readcore

   end module mod_coregeom
