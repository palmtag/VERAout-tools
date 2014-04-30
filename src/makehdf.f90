!=======================================================================
!
!  Program to create HDF test file with boron and eigenvalue at multiple statepoints
!
!  Copyright (c) 2014 Core Physics, Inc.
!
!  Distributed under the MIT license.  
!  See the LICENSE file in the main directory for details.
!
!  20140426 - modify to use lower case datasets
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
      integer     :: n

      integer     :: mapcore(3,3)
      data mapcore / 1, 2, 3, 4, 5, 6, 7, 8, 9 /

      real(8)     :: atemp(1), db
      real(8)     :: axial(10)

      real(8), allocatable :: xkeff(:)      ! eigenvalues
      real(8), allocatable :: xboron(:)     ! boron concentrations
      real(8), allocatable :: xexp(:)       ! exposure points (GWD/MT)

      character(len=10) :: cdate, ctime
      character(len=80) :: stringin      ! make these strings different lengths
      character(len=11) :: group_name

      integer     :: nstate
      integer     :: ierror


!--- read number of statepoints and create fake data

      write (*,*) 'How many statepoints?'
      read (*,*) nstate

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


!--- allocate and fill dummy arrays

!--- call wrapper routine to write integer 0D array

      idims(:)=0       ! clear

      group_name='/STATE_0000'    ! number will be overwritten below

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

! version number

      idims(1)=1
      itemp(1)=0       ! hack: send as single element array in order to pass type checking

      dsetname = 'veraout_version'
      call hwrite_integer(file_id, dsetname, idims, itemp)

! symmetry

      itemp(1)=0       ! hack: send as single element array in order to pass type checking

      dsetname = 'core_sym'
      call hwrite_integer(file_id, dsetname, idims, itemp)

! version

      itemp(1)=1       ! hack: send as single element array in order to pass type checking

      dsetname = 'version'
      call hwrite_integer(file_id, dsetname, idims, itemp)

! core map

      idims(:)=0     ! clear
      idims(1)=3
      idims(2)=3

      dsetname = 'core_map'
      call hwrite_integer(file_id, dsetname, idims, mapcore)

! axial mesh

      do n=1, 10
        axial(n)=dble(n)*10.0d0
      enddo

      idims(:)=0     ! clear
      idims(1)=10

      dsetname = 'axial_mesh'
      call hwrite_double(file_id, dsetname, idims, axial)


! title string

      dsetname='title'
      stringin='Sample HDF output file without pin data'
      call hwrite_string(file_id, dsetname, stringin)


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
! BUILD_DATE string string
! RUN_DATE   string string Required 
! RUN_TIME   string string Required
! PIN_POWERS   4 double
! PIN_VOLUMES  4 double
! PIN_EXPOSURE 4 double

      do n=1, nstate

! create group for statepoint

        write (group_name(8:11),'(i4.4)') n
        write (*,*) 'create statepoint ', n, group_name

        call h5gcreate_f (file_id, group_name, group_id, ierror)
        if (ierror.ne.0) stop 'error creating group'

! boron

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=xboron(n)   ! hack: send as single element array in order to pass type checking

        dsetname = group_name//'/boron'
        call hwrite_double(file_id, dsetname, idims, atemp)

! eigenvalue

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=xkeff(n)   ! hack: send as single element array in order to pass type checking

        dsetname = group_name//'/keff'
        call hwrite_double(file_id, dsetname, idims, atemp)

! exposure

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=xexp(n)   ! hack: send as single element array in order to pass type checking

        dsetname = group_name//'/EXPOSURE'
        call hwrite_double(file_id, dsetname, idims, atemp)

! power

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=100.0d0   ! hack: send as single element array in order to pass type checking

        dsetname = group_name//'/power'
        call hwrite_double(file_id, dsetname, idims, atemp)

! flow

        idims(:)=0     ! clear
        idims(1)=1
        atemp(1)=100.0d0   ! hack: send as single element array in order to pass type checking

        dsetname = group_name//'/flow'
        call hwrite_double(file_id, dsetname, idims, atemp)

!--- write strings

        dsetname=group_name//'/build_date'
        stringin= '2010-01-01'
        call hwrite_string(file_id, dsetname, stringin)

       call date_and_time (cdate, ctime)   ! , zone, values)

        dsetname=group_name//'/run_date'
        call hwrite_string(file_id, dsetname, cdate)

        dsetname=group_name//'/run_time'
        call hwrite_string(file_id, dsetname, ctime)

!---------------

! close group

        call h5gclose_f (group_id, ierror)
        if (ierror.ne.0) stop 'error closing group'

      enddo    ! end loop over statepoints

!-------------------------------
!   finished
!-------------------------------

!--- close file.

      call h5fclose_f(file_id, ierror)

      write (*,'(2a)') ' finished writing file: ', trim(filename)

      deallocate (xkeff)
      deallocate (xboron)
      deallocate (xexp)

!--- Close fortran interface.

      call h5close_f(ierror)

      end

