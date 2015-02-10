!=======================================================================
!
!  Module to include all subroutines for reading HDF5 file
!
!  Copyright (c) 2014 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!  Main directory: https://github.com/palmtag/VERAout-tools/
!
!=======================================================================
!
!  Routine to return type and dimensions:
!
!     h5info(file_id, dataset, itype, ndim, idim)
!       ** currently returns itype=0 always
!       ** currently does not work for character arrays
!
!  Routines to read:
!
!     hdf5_read_double  (generic subroutine)
!        read_double   (file_id, dataset, xvar)
!        read_double1d (file_id, dataset, max_idm, idm, xvar)
!        read_double2d (file_id, dataset, idm1, idm2, xvar)
!        read_double3d (file_id, dataset, idm1, idm2, idm3, xvar)
!        read_double4d (file_id, dataset, idm1, idm2, idm3, idm4, xvar)
!     read_string   (file_id, dataset, stringout)
!     read_string1d (file_id, dataset, stringout, maxsize, nsize)
!
!     hdf5_read_integer  (generic subroutine)
!        read_integer  (file_id, dataset, ivar)
!        read_integer1d(file_id, dataset, max_idm, idm, ivar)
!        read_integer2d(file_id, dataset, idm, jdm, ivar)
!
!  Routines to write:
!
!     hwrite_integer(file_id, dataset, idims, xvar)
!     hwrite_real   (file_id, dataset, idims, xvar)
!     hwrite_double (file_id, dataset, idims, xvar)
!     hwrite_string (file_id, dataset, strbuf)
!     hwrite_stringx(file_id, dataset, idims, xvar) **** fix: arrays of strings
!
!
!   Prefix  H5A Attributes
!     H5D   Datasets
!     H5E   Error reports
!     H5F   Files
!     H5G   Groups
!     H5I   Identifiers
!     H5L   Links
!     H5O   Objects
!     H5P   Property lists
!     H5R   References
!     H5S   Dataspaces
!     H5T   Datatypes
!     H5Z   Filters
!
!  Note: hid_t is 4 on the linux boxes I tested
!
!=======================================================================

      module mod_hdftools
      use     hdf5
      implicit none

      logical :: ifdebug=.false.    ! common debug flag for all subroutines

      private :: ifdebug

      private :: read_integer
      private :: read_integer1d
      private :: read_integer2d

      private :: read_double
      private :: read_double1d
      private :: read_double2d
      private :: read_double3d
      private :: read_double4d
      private :: h5_fatal

!--- define generic interface to read datasets

      interface hdf5_read_integer
        module procedure read_integer, read_integer1d, read_integer2d
      end interface hdf5_read_integer

      interface hdf5_read_double
        module procedure read_double, read_double1d, read_double2d, read_double3d, read_double4d
      end interface hdf5_read_double

! *** TO-DO: Create generic interface to write data

      contains

!=======================================================================
!
!  Subroutine to read 1D array of integers from HDF file
!
!  integer array has already been allocated, but can be larger than array to be read
!  (this is different than read_integer2d)
!
!=======================================================================
      subroutine read_integer1d(file_id, dataset, max_idm, idm, ivar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: max_idm
      integer,          intent(out):: idm
      integer,          intent(out):: ivar(max_idm)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: k

      integer, parameter :: max_dimen=2     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      idm=0    ! init
      ivar(:)=0

!--- open a dataset and dataspace

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

      call h5dget_space_f(dset_id, ispace_id, ierror)              !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f', 'ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

!  allow size 0 for SCALAR data

      if (maxndim.ne.1 .and. maxndim.ne.0) then
        write (*,*) 'invalid number of dimensions in read_integer1d'
        write (*,*) '  expecting 0 or 1, found ', maxndim
        call h5_fatal('read_integer1d', 'invalid number of dimensions')
      endif

      idm=int(h_dims(1))
      if (ifdebug) then
        write (*,*) 'idm ', idm
      endif
      if (idm.gt.max_idm) then
        write (*,*) 'max_idm = ', max_idm
        write (*,*) 'idm     = ', idm
        write (*,*) 'ERROR: max_idm exceeded'
        idm=0
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ivar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('h5dread_f','error reading data set '//trim(dataset))
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) then
        write (*,'(/,a)') ' 1D Integer Array'
        do k=1, idm
          write (*,'(2x,i4,20i12)') k, ivar(k)
        enddo
      endif

      return
      end subroutine read_integer1d

!=======================================================================
!
!  Subroutine to read 2D array of integers from HDF file
!
!  integer array has already been allocated and exact size must be given
!  (hint: find size from call to h5info and allocate array accordingly)
!
!  output:  ivar(idm,jdm)
!
!=======================================================================
      subroutine read_integer2d(file_id, dataset, idm, jdm, ivar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm, jdm
      integer,          intent(out):: ivar(idm,jdm)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: i, j

      integer, parameter :: max_dimen=2     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      ivar(:,:)=0

!--- open a dataset and dataspace

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

      call h5dget_space_f(dset_id, ispace_id, ierror)              !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

      if (maxndim.gt.max_dimen) call h5_fatal('read_integer2d','maximum dimensions exceeded')
      if (maxndim.ne.2) then
        write (*,*) 'invalid number of dimensions in read_integer2d'
        write (*,*) '  expecting 2, found ', maxndim
        call h5_fatal('read_integer2d','invalid number of dimensions')
      endif

      if (h_dims(1).ne.idm) call h5_fatal('read_integer2d','idm error')
      if (h_dims(2).ne.jdm) call h5_fatal('read_integer2d','jdm error')

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ivar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_integer2d','error reading data set '//trim(dataset))
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) then
        write (*,'(/,a)') ' 2D integer array:'
        do j=1, idm
          write (*,'(2x,i4,20i12)') j, (ivar(i,j),i=1,idm)
        enddo
      endif

      return
      end subroutine read_integer2d

!=======================================================================
!
!  Subroutine to read integer scalar from HDF5 file
!
!  This routine could be generalized further if we passed data type
!
!=======================================================================
      subroutine read_integer(file_id, dataset, ivar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(out):: ivar

!--- local

      integer :: idm_max
      integer :: idm
      integer :: iarray_temp(1)

      ivar=-1         ! init

!  integer scalars are actually 1D arrays with a single element

      idm_max=1
      idm=1

      call read_integer1d (file_id, dataset, idm_max, idm, iarray_temp)

      ivar=iarray_temp(1)

      return
      end subroutine read_integer

!=======================================================================
!
!  Subroutine to read double precision scalar from HDF5 file
!
!  This routine could be generalized further if we passed data type
!
!=======================================================================
      subroutine read_double(file_id, dataset, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      real(8),          intent(out):: xvar

!--- local

      integer :: idm_max
      integer :: idm
      real(8) :: array_temp(1)

      xvar=-1.0d0    ! init

!--- Read double precision scalar

!  double precisions scalars are actually 1D arrays with a single element

      idm_max=1
      idm=1

      call read_double1d (file_id, dataset, idm_max, idm, array_temp)

      xvar=array_temp(1)

      return
      end subroutine read_double
!=======================================================================
!
!  Subroutine to read 1D array of double from HDF file
!
!  double array has already been allocated, but can be larger than array to be read
!  (this is different than 2D read)
!
!=======================================================================
      subroutine read_double1d(file_id, dataset, max_idm, idm, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: max_idm
      integer,          intent(out):: idm
      real(8),          intent(out):: xvar(max_idm)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: k

      integer, parameter :: max_dimen=2     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      idm=0    ! init
      xvar(:)=0.0d0

!--- open a dataset and dataspace

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

      call h5dget_space_f(dset_id, ispace_id, ierror)              !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

!  allow size 0 for SCALAR data

      if (maxndim.gt.max_dimen) call h5_fatal('read_double1d','maximum dimensions exceeded')
      if (maxndim.ne.1 .and. maxndim.ne.0) then
        write (*,*) 'invalid number of dimensions in read_double1d'
        write (*,*) '  expecting 0 or 1, found ', maxndim
        call h5_fatal('read_double1d','invalid number of dimensions')
      endif

      idm =int(h_dims(1))
      if (ifdebug) then
        write (*,*) 'idm ', idm
      endif
      if (idm.gt.max_idm) then
        write (*,*) 'max_idm = ', max_idm
        write (*,*) 'idm     = ', idm
        write (*,*) 'ERROR: max_idm exceeded'
        idm=0
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xvar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal ('read_double1d', 'error reading data set '//trim(dataset))
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) then
        write (*,'(/,a)') ' 1D Double array:'
        do k=1, idm
          write (*,'(2x,i4,20f12.4)') k, xvar(k)
        enddo
      endif

      return
      end subroutine read_double1d
!=======================================================================
!
!  Subroutine to read 2D array of double from HDF file
!
!  double array has already been allocated and the dimensions
!  must match what is on the file exactly
!  (this is different than 1D read)
!
!=======================================================================
      subroutine read_double2d(file_id, dataset, idm1, idm2, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm1, idm2
      real(8),          intent(out):: xvar(idm1,idm2)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: i, j

      integer, parameter :: max_dimen=2     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      xvar(:,:)=0.0d0

!--- open a dataset and dataspace

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

      call h5dget_space_f(dset_id, ispace_id, ierror)              !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

      if (maxndim.gt.max_dimen) call h5_fatal('read_double2d','maximum dimensions exceeded')
      if (maxndim.ne.2) then
        write (*,*) 'invalid number of dimensions in read_double2d'
        write (*,*) '  expecting 2, found ', maxndim
        call h5_fatal('read_double2d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)).or. &
          idm2.ne.int(h_dims(2))) then
        write (*,*) 'input idm1 idm2 = ', idm1, idm2
        write (*,*) 'output hdim   = ', h_dims(1:2)
        write (*,*) 'ERROR: array bounds do not match'
        call h5_fatal('read_double2d','array bounds do not match')
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xvar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_double2d','error reading data set '//trim(dataset))
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) then
        write (*,'(/,a)') ' 2D Double array:'
        do j=1, idm2
          write (*,'(2x,i4,20f12.4)') j, (xvar(i,j),i=1,idm1)
        enddo
      endif

      return
      end subroutine read_double2d
!=======================================================================
!
!  Subroutine to read 3D array of double from HDF file
!
!  double array has already been allocated and the dimensions
!  must match what is on the file exactly
!  (this is different than 1D read)
!
!=======================================================================
      subroutine read_double3d(file_id, dataset, idm1, idm2, idm3, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm1, idm2, idm3
      real(8),          intent(out):: xvar(idm1,idm2,idm3)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim

      integer, parameter :: max_dimen=4     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      xvar(:,:,:)=0.0d0

!--- open a dataset and dataspace

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

      call h5dget_space_f(dset_id, ispace_id, ierror)              !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

      if (maxndim.gt.max_dimen) call h5_fatal('read_double3d','maximum dimensions exceeded')
      if (maxndim.ne.3) then
        write (*,*) 'invalid number of dimensions in read_double3d'
        write (*,*) '  expecting 3, found ', maxndim
        call h5_fatal('read_double3d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)) .or. &
          idm2.ne.int(h_dims(2)) .or. &
          idm3.ne.int(h_dims(3))) then
        write (*,*) 'input idm = ', idm1, idm2, idm3
        write (*,*) 'output hdim   = ', h_dims(1:3)
        write (*,*) 'ERROR: array bounds do not match'
        call h5_fatal('read_double3d','array bounds do not match')
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xvar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_double3d','error reading dataset')
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

!!*** too big to print
!!    if (ifdebug) then
!!      write (*,'(/,a)') ' 3D Double array:'
!!      do k=1, idm3
!!        do j=1, idm2
!!          write (*,'(2x,2i4,20f12.4)') k, j, (xvar(i,j,k),i=1,idm)
!!        enddo
!!      enddo
!!    endif

      return
      end subroutine read_double3d

!=======================================================================
!
!  Subroutine to read 4D array of double from HDF file
!
!  double array has already been allocated and the dimensions
!  must match what is on the file exactly
!  (this is different than 1D read)
!
!=======================================================================
      subroutine read_double4d(file_id, dataset, idm1, idm2, idm3, idm4, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm1, idm2, idm3, idm4
      real(8),          intent(out):: xvar(idm1,idm2,idm3,idm4)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim

      integer, parameter :: max_dimen=4     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      xvar(:,:,:,:)=0.0d0

!--- open a dataset and dataspace

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

      call h5dget_space_f(dset_id, ispace_id, ierror)              !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

      if (maxndim.gt.max_dimen) call h5_fatal('read_double4d','maximum dimensions exceeded')
      if (maxndim.ne.4) then
        write (*,*) 'invalid number of dimensions in read_double4d'
        write (*,*) '  expecting 4, found ', maxndim
        call h5_fatal('read_double4d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)) .or. &
          idm2.ne.int(h_dims(2)) .or. &
          idm3.ne.int(h_dims(3)) .or. &
          idm4.ne.int(h_dims(4))) then
        write (*,*) 'input idm = ', idm1, idm2, idm3, idm4
        write (*,*) 'output hdim   = ', h_dims(1:4)
        write (*,*) 'ERROR: array bounds do not match'
        call h5_fatal('read_double4d','array bounds do not match')
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xvar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_double4d','error reading dataset')
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

!!*** too big to print
!!    if (ifdebug) then
!!      write (*,'(/,a)') ' 4D Double array:'
!!      do j=1, idm2
!!        write (*,'(2x,i4,20f12.4)') j, (xvar(i,j),i=1,idm)
!!      enddo
!!    endif

      return
      end subroutine read_double4d
!=======================================================================
!
!  Subroutine to read character string from HDF5 file
!
!=======================================================================
      subroutine read_string(file_id, dataset, stringout)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      character(len=*), intent(out):: stringout

!--- local

      integer(hid_t)     :: dset_id
      integer(hid_t)     :: type_id
      integer            :: ierror

      integer, parameter :: max_dimen=2     ! maximum number of dimensions allowed

      integer(hsize_t)   :: h_dims(max_dimen)
      integer(size_t)    :: strlen
      integer(size_t)    :: iht

      ierror=0
      stringout=' '    ! initialize

!--- read string

!--- open a dataset

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

!--- determine data type (scalars do not use dataspace)

      call h5dget_type_f(dset_id, type_id, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_type_f ', ierror
        call h5_fatal('h5dget_type_f','ierror')
      endif

!  note that type_id does not match H5T_NATIVE_CHARACTER or H5T_STRING

      call h5tget_size_f(type_id, strlen, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_size_f ', ierror
        call h5_fatal('h5dget_size_f','ierror')
      endif

      if (strlen.gt.len(stringout)) then
        write (*,*) 'file  strlen   ', strlen
        write (*,*) 'input strlenin ', len(stringout)
        write (*,*) 'error: insufficient string length defined for read'
        call h5_fatal('read_string','insufficient string length defined for read')
      endif

! Read data from dataset

      call h5dread_f(dset_id, type_id, stringout, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dread_f reading data set ',trim(dataset)
        call h5_fatal('read_string','error reading dataset')
      endif

      do iht=1, strlen    ! search for null
        if (ichar(stringout(iht:iht)).eq.0) then
          stringout(iht:iht)=' '  ! remove null
        endif
      enddo

! Close datatype and dataset

      call h5tclose_f(type_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) write (*,'(1x,3a)') 'debug: read_string >',trim(stringout),'<'

      return
      end subroutine read_string

!=======================================================================
!
!  Subroutine to read 1D array of character string from HDF5 file
!
!=======================================================================
      subroutine read_string1d(file_id, dataset, stringout, maxsize, nsize)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: maxsize   ! maximum size of array
      integer,          intent(out):: nsize     ! returned size of array
      character(len=*), intent(out):: stringout(maxsize)

!--- local

      integer(hid_t)     :: dset_id
      integer(hid_t)     :: type_id
      integer(hid_t)     :: ispace_id
      integer            :: ierror, i, j, ilen
      integer            :: maxndim

      integer, parameter :: max_dimen=2     ! maximum number of dimensions allowed

      integer(hsize_t)   :: h_dims(max_dimen)
      integer(hsize_t)   :: h_maxdims(max_dimen)

      character(len=1), allocatable :: htemp(:,:)

      ierror=0
      nsize=0
      stringout(:)=' '    ! initialize

!--- read string

!--- open a dataset

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

!--- determine data type

      call h5dget_type_f(dset_id, type_id, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_type_f ', ierror
        call h5_fatal('h5dget_type_f','ierror')
      endif

!  note that type_id does not match H5T_NATIVE_CHARACTER or H5T_STRING

! open dataspace

      call h5dget_space_f(dset_id, ispace_id, ierror)       !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

      if (maxndim.gt.max_dimen) call h5_fatal('read_string1d','maximum dimensions exceeded')
      if (maxndim.ne.2) then
         write (*,*) 'maxndim  = ', maxndim
         write (*,*) 'expecting 2'
         write (*,*) 'max_dimen= ', max_dimen
         call h5_fatal('read_string1d','invalid number of dimensions')
      endif

!--- Note that the length of the character strings returned from the read
!--- are not the same as the lengths passed into the subroutine
!--- To fix, read data into a temporary character array then fill
!--- the actual character arrays

      ilen =int(h_dims(1))   ! temporary array length
      nsize=int(h_dims(2))
      allocate (htemp(ilen,nsize))

! Read data from dataset - read temporary character array

      call h5dread_f(dset_id, type_id, htemp, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dread_f reading data set ',trim(dataset)
        call h5_fatal('read_string1d','error reading dataset')
      endif

! Close datatype and dataset

      call h5tclose_f(type_id, ierror)
      call h5dclose_f(dset_id, ierror)

! Transfer temporary character arrays to actual character arrays

      if (ilen.gt.len(stringout)) then         ! protect bounds
        write (*,*) 'warning: truncating string length in read_string1d'
        ilen=len(stringout)
      endif
      if (nsize.gt.maxsize) then               ! protect bounds
        write (*,*) 'warning: truncating number of data entries in read_string1d'
        nsize=maxsize
      endif

      stringout(:)=' '
      do j=1, nsize
        do i=1, ilen
          if (ichar(htemp(i,j)).eq.0) htemp(i,j)=' '   ! remove null
          stringout(j)(i:i)=htemp(i,j)
        enddo
      enddo

      deallocate (htemp)

      if (ifdebug) then
        write (*,'(1x,2a)') 'dataset: ', trim(dataset)
        do i=1, nsize
          write (*,'(i6,3a)') i, '>',trim(stringout(i)),'<'
        enddo
      endif

      return
      end subroutine read_string1d

!=======================================================================
!
!  Wrapper Subroutine to write a real array to an HDF5 file
!
!  file_id - file id returned by HDF5 open statement
!  dataset - data set name
!  idims   - array that specifies the size of each dimension
!       set any unused dimensions to zero
!       examples:
!         idims= 100,  0,  0, 0, 0, 0, 0, 0, 0, 0   for 1D array
!         idims= 100, 20,  0, 0, 0, 0, 0, 0, 0, 0   for 2D array
!         idims= 100, 20, 16, 0, 0, 0, 0, 0, 0, 0   for 3D array, etc.
!  xvar    - data array to be written
!
!=======================================================================
      subroutine hwrite_real(file_id, dataset, idims, xvar)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: idims(10)    ! Dataset dimensions
      real,             intent(in) :: xvar(*)      ! multi-dimensional array, treated as 1D

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(hid_t) :: dtype_id      ! data type to be written

      integer  :: i
      integer  :: irank               ! Dataset rank
      integer  :: ierror              ! Error flag
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

!--- calculate rank from idims array

      dtype_id=H5T_NATIVE_REAL

      irank=0
      do i=1, 10
        idimht(i)=idims(i)    ! convert from normal integer to hsize_t
        if (idimht(i).gt.0) irank=i
      enddo

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset
      write (*,30) 'rank ', irank
      write (*,30) 'dims ', idimht(1:irank)
   30 format (3x,'creating dataspace with ',a, 10i6)

! Create the dataspace (information about array)

      call h5screate_simple_f(irank, idimht, dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5screate_simple_f','ierror')

! Create the dataset with default properties.

      call h5dcreate_f(file_id, dataset, dtype_id, dspace_id, dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5dcreate_f','ierror')

! Write dataset

      call h5dwrite_f(dset_id, dtype_id, xvar, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('h5dwrite_f','ierror')

! End access to the dataset and release resources used by it. (release dataset id)

      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5dclose_f','ierror')

! Terminate access to the data space.  (release dataspace id)

      call h5sclose_f(dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5sclose_f','ierror')

      return
      end subroutine hwrite_real
!=======================================================================
!
!  Wrapper Subroutine to write a double array to an HDF5 file
!
!  file_id - file id returned by HDF5 open statement
!  dataset - data set name
!  idims   - array that specifies the size of each dimension
!       set any unused dimensions to zero
!       examples:
!         idims= 100,  0,  0, 0, 0, 0, 0, 0, 0, 0   for 1D array
!         idims= 100, 20,  0, 0, 0, 0, 0, 0, 0, 0   for 2D array
!         idims= 100, 20, 16, 0, 0, 0, 0, 0, 0, 0   for 3D array, etc.
!  xvar    - data array to be written
!
!=======================================================================
      subroutine hwrite_double(file_id, dataset, idims, xvar)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: idims(10)    ! Dataset dimensions
      real(8),          intent(in) :: xvar(*)      ! multi-dimensional array, treated as 1D

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(hid_t) :: dtype_id      ! data type to be written

      integer  :: i
      integer  :: irank               ! Dataset rank
      integer  :: ierror              ! Error flag
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

!--- calculate rank from idims array

      dtype_id=H5T_NATIVE_DOUBLE

      do i=1, 10
        idimht(i)=idims(i)    ! convert from normal integer to hsize_t
        if (idimht(i).gt.0) irank=i
      enddo

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset
      write (*,30) 'rank ', irank
      write (*,30) 'dims ', idimht(1:irank)
   30 format (3x,'creating dataspace with ',a, 10i6)

! Create the dataspace (information about array)

!*** to create a true Scalar value (not a 1D array with one element), 
!*** you could use "h5screate_f(H5S_SCALAR_F, dspace_id, ierror)   (not verified)

      call h5screate_simple_f(irank, idimht, dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5screate_simple','ierror')

! Create the dataset with default properties.

      call h5dcreate_f(file_id, dataset, dtype_id, dspace_id, dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5dcreate_f','ierror')

! Write dataset

      call h5dwrite_f(dset_id, dtype_id, xvar, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('h5dwrite_f','ierror')

! End access to the dataset and release resources used by it. (release dataset id)

      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5dclose_f','ierror')

! Terminate access to the data space.  (release dataspace id)

      call h5sclose_f(dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5sclose_f','ierror')


      return
      end subroutine hwrite_double

!=======================================================================
!
!  Wrapper Subroutine to write a integer array to an HDF5 file
!
!  file_id - file id returned by HDF5 open statement
!  dataset - data set name
!  idims   - array that specifies the size of each dimension
!       set any unused dimensions to zero
!       examples:
!         idims= 100,  0,  0, 0, 0, 0, 0, 0, 0, 0   for 1D array
!         idims= 100, 20,  0, 0, 0, 0, 0, 0, 0, 0   for 2D array
!         idims= 100, 20, 16, 0, 0, 0, 0, 0, 0, 0   for 3D array, etc.
!  xvar    - data array to be written
!
!=======================================================================
      subroutine hwrite_integer(file_id, dataset, idims, xvar)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: idims(10)    ! Dataset dimensions
      integer,          intent(in) :: xvar(*)      ! multi-dimensional array, treated as 1D

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(hid_t) :: dtype_id      ! data type to be written

      integer  :: i
      integer  :: irank               ! Dataset rank
      integer  :: ierror              ! Error flag
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

!--- calculate rank from idims array

      dtype_id=H5T_NATIVE_INTEGER

      irank=0
      do i=1, 10
        idimht(i)=idims(i)    ! convert from normal integer to hsize_t
        if (idimht(i).gt.0) irank=i
      enddo

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset
      write (*,30) 'rank ', irank
      write (*,30) 'dims ', idimht(1:irank)
   30 format (3x,'creating dataspace with ',a, 10i6)

! Create the dataspace (information about array)

!*** to create a true Scalar value (not a 1D array with one element), 
!*** you could use "h5screate_f(H5S_SCALAR_F, dspace_id, ierror)   (not verified)

      call h5screate_simple_f(irank, idimht, dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5screate_simple','ierror')

! Create the dataset with default properties.

      call h5dcreate_f(file_id, dataset, dtype_id, dspace_id, dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dcreate_f','ierror')

! Write dataset

      call h5dwrite_f(dset_id, dtype_id, xvar, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dwrite_f','ierror')

! End access to the dataset and release resources used by it. (release dataset id)

      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dclose_f','ierror')

! Terminate access to the data space.  (release dataspace id)

      call h5sclose_f(dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5sclose_f','ierror')

      return
      end subroutine hwrite_integer

!=======================================================================
!
!  Wrapper Subroutine to write a character string to an HDF5 file
!
!  file_id - file id returned by HDF5 open statement
!  dataset - data set name
!  strbuf  - character string to be written
!
!=======================================================================
      subroutine hwrite_string(file_id, dataset, strbuf)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      character(len=*), intent(in) :: strbuf       ! character string

!--- local

      integer(hid_t)   :: dset_id       ! Dataset identifier
      integer(hid_t)   :: dspace_id     ! Dataspace identifier
      integer(hid_t)   :: type_id       ! Datatype identifier
      integer(hsize_t) :: dimstr(1)
      integer(size_t)  :: lenstr        ! fixed - was integer*8 - fixed again, was integer

      integer ::   ierror        ! Error flag

      lenstr=len_trim(strbuf)   ! len(strbuf)  write actual length, not character length - won't work for arrays
      dimstr(1)=1

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset

      call h5tcopy_f (H5T_NATIVE_CHARACTER, type_id, ierror)
!!    print *, 'h5tcopy_f returns type_id', type_id

      call h5tset_size_f (type_id, lenstr, ierror)
!!    print *, 'h5tset_size_f returns ierror', ierror

!--- Create the dataspace.

      call h5screate_f (H5S_SCALAR_F, dspace_id, ierror)
!!    print *, 'h5screate_f returns dspace_id', dspace_id

!--- Create the dataset with default properties.

      call h5dcreate_f(file_id, dataset, type_id, dspace_id, dset_id, ierror)
!!    print *, 'h5dcreate_f returns', dset_id

      call h5dwrite_f (dset_id, type_id, strbuf, dimstr, ierror)
!!    print *, 'h5dwrite_f returns', ierror

      call h5tclose_f(type_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5tclose_f','ierror')

      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dclose_f','ierror')

      call h5sclose_f(dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5sclose_f','ierror')

      return
      end subroutine hwrite_string
!=======================================================================
!
!  Wrapper Subroutine to write a character string to an HDF5 file
!
!  file_id - file id returned by HDF5 open statement
!  dataset - data set name
!
! ***** fix: this is not working, how do you write an array of strings? *******
!
!=======================================================================
      subroutine hwrite_stringx(file_id, dataset, idims, namex)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: idims(10)    ! Dataset dimensions
      character(len=*), intent(in) :: namex(*)     ! multi-dimensional array, treated as 1D

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
!att  integer(hid_t) :: aspace_id     ! attribute identifier
!att  integer(hid_t) :: attr_id       ! attribute id
      integer(size_t):: attrlen       ! length of string
      integer(hid_t) :: type_id       ! attribute type to be written

      integer  :: i
      integer  :: ierror              ! Error flag
      integer  :: irank               ! dataset rank
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

! Create scalar dataspace (information about array)

      attrlen=len(namex)
      irank=0
      do i=1, 10
        idimht(i)=idims(i)    ! convert from normal integer to hsize_t
        if (idimht(i).gt.0) irank=i
      enddo

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset
      write (*,*)         ' debug: rank           : ', irank
      write (*,*)         ' debug: size           : ', idimht(1:irank)
      write (*,*)         ' debug: len            : ', attrlen
!att  call h5screate_simple_f(irank, idimht, aspace_id, ierror)
!att  if (ierror.ne.0) call h5_fatal('ierror: h5screate_simple','ierror')

      write (*,*) '*** WARNING: the ability to write arrays of strings is not working ***'   ! fix
      write (*,*) '*** WARNING: This needs to be fixed                                ***'   ! fix

!-------------------------------------------------------
!  For an example of writing a 1D array of strings,
!  see example input "attrexample.f90"
!-------------------------------------------------------

! Create datatype for the attribute.

      call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierror)
      call h5tset_size_f(type_id, attrlen, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5tset_size_f','ierror')

! Create the dataspace

      call h5screate_f (H5S_SCALAR_F, dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5screate_f','ierror')

! Create the dataset with default properties.

      call h5dcreate_f(file_id, dataset, type_id, dspace_id, dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dcreate_f','ierror')

! Create dataset attribute.

      call h5dwrite_f (dset_id, type_id, namex, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dwrite_f','ierror')

!att  call h5acreate_f(dset_id, dataset, type_id, aspace_id, attr_id, ierror)
!att  if (ierror.ne.0) call h5_fatal('ierror: h5acreate_f','ierror')

! Write the attribute data.

!att  call h5awrite_f(attr_id, type_id, namex, idimht, ierror)

! Close the attribute.

!att  call h5aclose_f(attr_id, ierror)

! Close the dataset

      call h5tclose_f(type_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5tclose_f','ierror')
      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dclose_f','ierror')
      call h5sclose_f(dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5sclose_f','ierror')

      return
      end subroutine hwrite_stringx
!=======================================================================
!
!  Subroutine to return size of array from HDF file
!
!  input:
!     file_id     file_id of an already opened file
!     dataset     dataset name to look for
!
!  output:
!     itype       dataset type - returns -99 if dataset not found
!                   **** returns 0 otherwise  **** (needs to be fixed)
!     ndim        number of dimensions - returns 0 for scalar
!     idim(10)    integer array of dimension sizes
!
!=======================================================================
      subroutine h5info(file_id, dataset, itype, ndim, idim)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset

      integer, intent(out) :: itype
      integer, intent(out) :: ndim
      integer, intent(out) :: idim(10)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer(hid_t) :: type_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: i
      logical        :: ifxst

      integer, parameter :: max_dimen=10    ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

!---

      itype=0
      ndim=0
      idim(:)=0
      ierror=0

      write (*,'(/,1x,2a)') 'checking dataset: ', trim(dataset)

!--- check if dataset exists

      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (.not.ifxst) then
        itype=-99
        return
      endif

!--- open a dataset and dataspace

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror<0) then
        write (*,*) 'error: data set ',trim(dataset), ' could not be opened'
        write (*,*) 'error code ', ierror
        return
      endif

!--- determine data type (scalars do not use dataspace)

      call h5dget_type_f(dset_id, type_id, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_type_f ', ierror
        call h5_fatal('h5dget_type_f', 'ierror')
      endif

!*** the following does not work
!     if     (type_id.eq.H5T_NATIVE_DOUBLE) then
!       itype=3
!     elseif (type_id.eq.H5T_NATIVE_REAL) then
!       itype=2
!     elseif (type_id.eq.H5T_NATIVE_INTEGER) then
!       itype=1
!     else
!       itype=99   ! unknown
!     endif

!--- open space

      call h5dget_space_f(dset_id, ispace_id, ierror)       !Assign a dataspace based on the dataset
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_space_f ', ierror
        call h5_fatal('error: h5dget_space_f','ierror')
      endif

! Determine number of discrete subunits for each dimensions
      h_dims(:)=0
      h_maxdims(:)=0
      call h5sget_simple_extent_dims_f(ispace_id, h_dims, h_maxdims, ierror)
      if (ierror.lt.0) then    ! returns rank in ierror
        write (*,*) 'error: h5sget_simple_extent_dims_f ', ierror
        call h5_fatal('error: h5sget_simple_extent_dims_f','ierror')
      endif
      maxndim=ierror

      if (maxndim.gt.max_dimen) call h5_fatal('h5info','maximum dimensions exceeded')

      if (ifdebug) then
        write (*,110) trim(dataset)
        write (*,115) 'number of dimensions: ', maxndim
        write (*,115) 'dimensions    : ', h_dims(1:maxndim)
        write (*,115) 'dimensions max: ', h_maxdims(1:maxndim)
      endif
  110 format (/,' dataset : ', a)
  115 format (1x,a,20i4)

      ndim=maxndim
      do i=1, ndim
        idim(i)=int(h_dims(i))   ! convert to standard integer
      enddo

!--- close datasets

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      return
      end subroutine h5info

!=======================================================================
!
!   Subroutine to call if fatal error
!
!   DO NOT USE STOP - Regression tests need an error code returned
!
!=======================================================================
      subroutine h5_fatal(subr,label)
      implicit none
      character(len=*) :: subr, label

      write (*,*) 'Fatal error in HDF subroutine ', trim(subr)
      write (*,*) 'Error: ', trim(label)

! use call to exit instead of stop to return an error code

      call exit(80)       ! return error code

      end subroutine h5_fatal
!=======================================================================

      end module mod_hdftools
