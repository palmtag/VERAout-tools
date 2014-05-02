!=======================================================================
!
!  Program to create VERAOUT test file with boron and eigenvalue at multiple statepoints
!
!  Copyright (c) 2014 Core Physics, Inc.
!
!  Distributed under the MIT license.  
!  See the LICENSE file in the main directory for details.
!
!  2014/04/26 - modify to use lower case datasets
!  2014/05/01 - update to match VERAOUT Rev. 1 
!
!=======================================================================
      program makehdf 

      use hdf5
      use mod_hdftools

      implicit none

      character(len=20) :: filename = 'output.h5'  ! File name
      character(len=40) :: dsetname                ! Dataset name

      integer(hid_t)   :: file_id       ! File identifier
      integer(hid_t)   :: group_id      ! Group identifier
      integer(hsize_t) :: idims(10)     ! Dataset dimensions (2D array)

      integer     :: itemp(10)
      integer     :: i, j, k, n, na

      integer :: nassm=9
      integer :: naxial=10
      integer :: npin=3
      real(8), allocatable :: pinpower(:,:,:,:)


      integer     :: mapcore(3,3)
      data mapcore / 1, 2, 3, 4, 5, 6, 7, 8, 9 /

      real(8)     :: atemp(1), db
      real(8)     :: xran   ! random number
      real(8)     :: sum

      real(8), allocatable :: axial(:)      ! axial mesh (naxial+1)
      real(8), allocatable :: xkeff(:)      ! eigenvalues
      real(8), allocatable :: xboron(:)     ! boron concentrations
      real(8), allocatable :: xexp(:)       ! exposure points (GWD/MT)

      character(len=10) :: cdate, ctime
      character(len=80) :: stringin      ! make these strings different lengths
      character(len=12) :: group_name

      integer     :: nstate
      integer     :: ierror


!--- read number of statepoints from input

      write (*,*) 'How many statepoints?'
      read (*,*) nstate

!--- allocate and fill dummy arrays

      allocate (xkeff(nstate))
      allocate (xboron(nstate))
      allocate (xexp(nstate))

      if (nstate.gt.1) then
        db=1000.0d0/dble(nstate-1)
      else
        db=0.0d0
      endif

      do n=1, nstate
        xexp(n)=dble(n-1)
        xkeff(n)=1.0d0 - (n-1)*1.0d-4
        xboron(n)=1000.0d0 - (n-1)*db
      enddo

      allocate (axial(naxial+1))
      allocate (pinpower(nassm, naxial, npin, npin))
      pinpower=0.0d0

      do n=1, naxial+1
        axial(n)=dble(n)*10.0d0
      enddo


!-------------------------------
!   Create HDF file
!-------------------------------

!--- Initialize fortran interface.

      call h5open_f(ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

      write (*,'(2a,i10)') 'creating file: ', trim(filename), file_id

!--- Create a new file using default properties.

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

!------------------------
!  Write Top-level data
!------------------------

! INPUT_PARAMETER_LIST group (various)
! VERAOUT_VERSION 0 integer Required
! CORE_MAP        2 integer Required
! CORE_SYM        0 integer Required
! TITLE           0 string Required

! version number - Version 1 corresponds to Rev. 1 specification

      idims(:)=0       ! clear
      idims(1)=1
      itemp(1)=1       ! hack: send as single element array in order to pass type checking

      dsetname = 'veraout_version'
      call hwrite_integer(file_id, dsetname, idims, itemp)

! title string

      dsetname='title'
      stringin='Sample HDF output file without pin data'
      call hwrite_string(file_id, dsetname, stringin)

!--- date stamps

      dsetname='build_date'
      stringin= '2010-01-01'
      call hwrite_string(file_id, dsetname, stringin)

      call date_and_time (cdate, ctime)   ! , zone, values)

      dsetname='run_date'
      call hwrite_string(file_id, dsetname, cdate)

      dsetname='run_time'
      call hwrite_string(file_id, dsetname, ctime)

!------------------------
!  Write CORE data
!------------------------

! create group for core data

      group_name='/CORE/'
      write (*,*) 'create group', group_name

      call h5gcreate_f (file_id, group_name, group_id, ierror)
      if (ierror.ne.0) stop 'error creating group'

! core map

      idims(:)=0     ! clear
      idims(1)=3
      idims(2)=3

      dsetname = trim(group_name)//'core_map'
      call hwrite_integer(file_id, dsetname, idims, mapcore)

! symmetry - 1=full core

      idims(:)=0     ! clear
      idims(1)=1
      itemp(1)=1     ! hack: send as single element array in order to pass type checking

      dsetname = trim(group_name)//'core_sym'
      call hwrite_integer(file_id, dsetname, idims, itemp)

! axial mesh

      idims(:)=0     ! clear
      idims(1)=naxial+1

      dsetname = trim(group_name)//'axial_mesh'
      call hwrite_double(file_id, dsetname, idims, axial)

! rated power and flow

      idims(:)=0     ! clear
      idims(1)=1
      atemp(1)=3400.0d0  ! hack: send as single element array in order to pass type checking

      dsetname = trim(group_name)//'rated_power'
      call hwrite_double(file_id, dsetname, idims, atemp)

      idims(:)=0     ! clear
      idims(1)=1
      atemp(1)=1800.0d0  ! hack: send as single element array in order to pass type checking

      dsetname = trim(group_name)//'rated_flow'
      call hwrite_double(file_id, dsetname, idims, atemp)

! assembly pitch

      idims(:)=0     ! clear
      idims(1)=1
      atemp(1)=21.54d0   ! hack: send as single element array in order to pass type checking

      dsetname = trim(group_name)//'apitch'
      call hwrite_double(file_id, dsetname, idims, atemp)

!  close CORE group

      call h5gclose_f (group_id, ierror)
      if (ierror.ne.0) stop 'error closing group'

!--------------------------
!  Write Statepoint data
!--------------------------

! AXIAL_MESH 1 double Required
! BORON      0 double
! KEFF       0 double
! POWER      0 double Required*
! FLOW       0 double
! TINLET     0 double
! EXPOSURE   0 double
! PIN_POWERS   4 double
! PIN_VOLUMES  4 double
! PIN_EXPOSURE 4 double

      group_name='/STATE_0000/'    ! number will be overwritten below

      do n=1, nstate

! create group for statepoint

        write (group_name(8:11),'(i4.4)') n
        write (*,*) 'create statepoint ', n, trim(group_name)

        call h5gcreate_f (file_id, group_name, group_id, ierror)
        if (ierror.ne.0) stop 'error creating group'

! boron

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=xboron(n)   ! hack: send as single element array in order to pass type checking

        dsetname = trim(group_name)//'boron'
        call hwrite_double(file_id, dsetname, idims, atemp)

! eigenvalue

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=xkeff(n)   ! hack: send as single element array in order to pass type checking

        dsetname = trim(group_name)//'keff'
        call hwrite_double(file_id, dsetname, idims, atemp)

! exposure

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=xexp(n)   ! hack: send as single element array in order to pass type checking

        dsetname = trim(group_name)//'exposure'
        call hwrite_double(file_id, dsetname, idims, atemp)

! power

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=100.0d0   ! hack: send as single element array in order to pass type checking

        dsetname = trim(group_name)//'power'
        call hwrite_double(file_id, dsetname, idims, atemp)

! flow

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=100.0d0   ! hack: send as single element array in order to pass type checking

        dsetname = trim(group_name)//'flow'
        call hwrite_double(file_id, dsetname, idims, atemp)


!--- pin power

        sum=0.0d0 
        do i=1, npin
          do j=1, npin
            if (i.eq.2 .and. j.eq.2) cycle   ! skip middle pin
            do k=1, naxial
              do na=1, nassm
                 call random_number(xran)
                 pinpower(na,k,j,i)= 1.0d0+xran
                 sum=sum+1.0d0+xran    ! assume uniform axial grid
              enddo
            enddo
          enddo
        enddo
        sum=nassm*naxial*(npin*npin-1)/sum
        write (*,*) 'debug: pin normalization=', sum

        do i=1, npin
          do j=1, npin
            if (i.eq.2 .and. j.eq.2) cycle   ! skip middle pin
            do k=1, naxial
              do na=1, nassm
                 pinpower(na,k,j,i)=pinpower(na,k,j,i)*sum
              enddo
            enddo
          enddo
        enddo

        idims(:)=0     ! clear
        idims(1)=nassm
        idims(2)=naxial
        idims(3)=npin
        idims(4)=npin

        dsetname = trim(group_name)//'pin_powers'
        call hwrite_double(file_id, dsetname, idims, pinpower)

!---------------

! close statepoint group

        call h5gclose_f (group_id, ierror)
        if (ierror.ne.0) stop 'error closing group'

      enddo    ! end loop over statepoints

!-------------------------------
!   finished
!-------------------------------

!--- close file

      call h5fclose_f(file_id, ierror)

      write (*,'(2a)') ' finished writing file: ', trim(filename)

      deallocate (xkeff)
      deallocate (xboron)
      deallocate (xexp)

!--- Close fortran interface.

      call h5close_f(ierror)

      end

