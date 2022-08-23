   module mod_hdftools
   use     hdf5
   implicit none
!=======================================================================
!
!  Module to include all wrapper subroutines for reading HDF5 file
!
!  Copyright (c) 2014-2017 Core Physics, Inc.
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
!        read_double1do(file_id, dataset, max_idm, idm, xvar)     ! exact size can be over allocated
!        read_double1d (file_id, dataset, idm1, xvar)
!        read_double2d (file_id, dataset, idm1, idm2, xvar)
!        read_double3d (file_id, dataset, idm1, idm2, idm3, xvar)
!        read_double4d (file_id, dataset, idm1, idm2, idm3, idm4, xvar)
!        read_double5d (file_id, dataset, idm1, idm2, idm3, idm4, idm5, xvar)
!     read_string   (file_id, dataset, stringout)
!     read_string1d (file_id, dataset, stringout, nsize)
!
!     hdf5_read_integer  (generic subroutine)
!        read_integer  (file_id, dataset, ivar)
!        read_integer1do(file_id, dataset, max_idm, idm, ivar)   ! exact size can be over allocated
!        read_integer1d(file_id, dataset, idm, ivar)
!        read_integer2d(file_id, dataset, idm, jdm, ivar)
!
!  Routines to write:
!
!     hwrite_integer(file_id, dataset, idims, xvar)
!     hwrite_integer_scalar(file_id, dataset, xvar)
!     hwrite_real   (file_id, dataset, idims, xvar)
!     hwrite_double (file_id, dataset, idims, xvar)
!     hwrite_double_scalar (file_id, dataset, xvar)
!     hwrite_string (file_id, dataset, strbuf)
!     hwrite_string1d(file_id, dataset, xvar, idims)  !  1D array of strings
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
!  Note: hid_t is equal to 4 on the linux boxes I tested
!
!
!=======================================================================
!
!  8/3/2016 - updated with optional argument to hwrite_double to allow user to update
!             a dataset instead of creating a new one
!             Need to port this ability to other write subroutines!!!
!
!  8/8/2017 - added hwrite_double_scalar
!           - added hwrite_integer_scalar
!            *** still need to write real scalars
!            *** and should write a generic interface
!
!  10/28/2017 - added read_double5d
!
!=======================================================================


      logical :: ifdebug=.false.    ! common debug flag for all subroutines

      private :: ifdebug

      private :: read_integer
      private :: read_integer1do
      private :: read_integer1d
      private :: read_integer2d

      private :: read_double
      private :: read_double1do
      private :: read_double1d
      private :: read_double2d
      private :: read_double3d
      private :: read_double4d
      private :: read_double5d
      private :: h5_fatal

!--- define generic interface to read datasets

      interface hdf5_read_integer
        module procedure read_integer, read_integer1do, read_integer1d, read_integer2d
      end interface hdf5_read_integer

      interface hdf5_read_double
        module procedure read_double, read_double1do, read_double1d, read_double2d, read_double3d, read_double4d, read_double5d
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
      subroutine read_integer1do(file_id, dataset, max_idm, idm, ivar)
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
        write (*,*) 'invalid number of dimensions in read_integer1do'
        write (*,*) '  expecting 0 or 1, found ', maxndim
        call h5_fatal('read_integer1do', 'invalid number of dimensions')
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
      end subroutine read_integer1do

!=======================================================================
!
!  Subroutine to read 1D array of integers from HDF file
!
!  integer array has already been allocated and exact size must be given
!  (hint: find size from call to h5info and allocate array accordingly)
!
!  output:  ivar(idm)
!
!=======================================================================
      subroutine read_integer1d(file_id, dataset, idm, ivar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm
      integer,          intent(out):: ivar(idm)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: i

      integer, parameter :: max_dimen=1     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_integer1d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_integer1d','invalid number of dimensions')
      endif

      if (h_dims(1).ne.idm) then
         write (0,*) '> number dimensions ', maxndim
         write (0,*) '> h_dims(1)         ', h_dims(1)
         write (0,*) '> idm               ', idm
         call h5_fatal('read_integer1d','idm error 1D')
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ivar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_integer1d','error reading data set '//trim(dataset))
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) then
        write (*,'(/,a)') ' 1D integer array:'
        do i=1, idm
          write (*,'(2x,i4,20i12)') i, ivar(i)
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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_integer2d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_integer2d','invalid number of dimensions')
      endif

      if (h_dims(1).ne.idm .or. h_dims(2).ne.jdm) then
        write (0,*) '> h_dims(1) ', h_dims(1)
        write (0,*) '> h_dims(2) ', h_dims(2)
        write (0,*) '> idm       ', idm
        write (0,*) '> jdm       ', jdm
        call h5_fatal('read_integer2d','idm-jdm error 2D')
      endif

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

!  integer scalars can be 1D arrays with a single element
!  The "1do" routine supports reading both ways

      idm_max=1
      idm=1

      call read_integer1do(file_id, dataset, idm_max, idm, iarray_temp)

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

!  double scalars can be 1D arrays with a single element
!  The "1do" routine supports reading both ways

      idm_max=1
      idm=1

      call read_double1do(file_id, dataset, idm_max, idm, array_temp)

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
      subroutine read_double1do(file_id, dataset, max_idm, idm, xvar)
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

      if (maxndim.ne.1 .and. maxndim.ne.0) then
        write (*,*) 'invalid number of dimensions in read_double1do'
        write (*,*) '  expecting 0 or 1, found ', maxndim
        call h5_fatal('read_double1do','invalid number of dimensions')
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
        call h5_fatal ('read_double1do', 'error reading data set '//trim(dataset))
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
      end subroutine read_double1do
!=======================================================================
!
!  Subroutine to read 2D array of double from HDF file
!
!  double array has already been allocated and the dimensions
!  must match what is on the file exactly
!  (this is different than 1D read)
!
!=======================================================================
      subroutine read_double1d(file_id, dataset, idm1, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm1
      real(8),          intent(out):: xvar(idm1)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim
      integer        :: i

      integer, parameter :: max_dimen=1     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_double1d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_double1d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1))) then
        write (*,*) 'input idm1  = ', idm1
        write (*,*) 'output hdim = ', h_dims(:)
        write (*,*) 'ERROR: array bounds do not match'
        call h5_fatal('read_double1d','array bounds do not match')
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xvar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_double1d','error reading data set '//trim(dataset))
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

      if (ifdebug) then
        write (*,'(/,a)') ' 1D Double array:'
        do i=1, idm1
          write (*,'(2x,i4,20f12.4)') i, xvar(i)
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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_double2d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_double2d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)).or. &
          idm2.ne.int(h_dims(2))) then
        write (*,*) 'input idm   = ', idm1, idm2
        write (*,*) 'output hdim = ', h_dims(:)
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

      integer, parameter :: max_dimen=3     ! maximum number of dimensions allowed

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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_double3d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_double3d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)) .or. &
          idm2.ne.int(h_dims(2)) .or. &
          idm3.ne.int(h_dims(3))) then
        write (*,*) 'input idm   = ', idm1, idm2, idm3
        write (*,*) 'output hdim = ', h_dims(:)
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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_double4d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_double4d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)) .or. &
          idm2.ne.int(h_dims(2)) .or. &
          idm3.ne.int(h_dims(3)) .or. &
          idm4.ne.int(h_dims(4))) then
        write (*,*) 'input idm   = ', idm1, idm2, idm3, idm4
        write (*,*) 'output hdim = ', h_dims(:)
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
!  Subroutine to read 5D array of double from HDF file
!
!  double array has already been allocated and the dimensions
!  must match what is on the file exactly
!  (this is different than 1D read)
!
!=======================================================================
      subroutine read_double5d(file_id, dataset, idm1, idm2, idm3, idm4, idm5, xvar)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: idm1, idm2, idm3, idm4, idm5
      real(8),          intent(out):: xvar(idm1,idm2,idm3,idm4, idm5)

!--- local

      integer(hid_t) :: dset_id
      integer(hid_t) :: ispace_id
      integer        :: ierror
      integer        :: maxndim

      integer, parameter :: max_dimen=5     ! maximum number of dimensions allowed

      integer(hsize_t) :: h_dims(max_dimen)
      integer(hsize_t) :: h_maxdims(max_dimen)

      ierror=0
      xvar(:,:,:,:,:)=0.0d0

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

      if (maxndim.ne.max_dimen) then
        write (*,*) 'invalid number of dimensions in read_double5d'
        write (*,*) '  expecting ', max_dimen
        write (*,*) '  found     ', maxndim
        call h5_fatal('read_double5d','invalid number of dimensions')
      endif

      if (idm1.ne.int(h_dims(1)) .or. &
          idm2.ne.int(h_dims(2)) .or. &
          idm3.ne.int(h_dims(3)) .or. &
          idm4.ne.int(h_dims(4)) .or. &
          idm5.ne.int(h_dims(5))) then
        write (*,*) 'input idm   = ', idm1, idm2, idm3, idm4, idm5
        write (*,*) 'output hdim = ', h_dims(:)
        write (*,*) 'ERROR: array bounds do not match'
        call h5_fatal('read_double5d','array bounds do not match')
        return
      endif

! Read data from dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xvar, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error reading data set ',trim(dataset)
        call h5_fatal('read_double5d','error reading dataset')
      endif

      call h5sclose_f(ispace_id, ierror)
      call h5dclose_f(dset_id, ierror)

!!*** too big to print
!!    if (ifdebug) then
!!      write (*,'(/,a)') ' 5D Double array:'
!!      do j=1, idm2
!!        write (*,'(2x,i4,20f12.4)') j, (xvar(i,j),i=1,idm)
!!      enddo
!!    endif

      return
      end subroutine read_double5d
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
      subroutine read_string1d(file_id, dataset, stringout, nsize)
      implicit none

      integer(hid_t),   intent(in) :: file_id
      character(len=*), intent(in) :: dataset
      integer,          intent(in) :: nsize     ! size of array
      character(len=*), intent(out):: stringout(nsize)

!--- local

      integer(hid_t)     :: dset_id
      integer(hid_t)     :: type_id
      integer(hid_t)     :: ispace_id
      integer            :: ierror, i, j
      integer            :: klen
      integer            :: maxndim

      integer, parameter :: max_dimen=1     ! maximum number of dimensions allowed

      integer(hsize_t)   :: h_dims(4)     ! (max_dimen)
      integer(hsize_t)   :: h_maxdims(4)  ! (max_dimen)

      character(len=1), allocatable :: htemp(:,:)

      ierror=0
      stringout(:)=' '    ! initialize

      klen =len(stringout)

!--- read string

!--- open a dataset

      write (*,'(1x,2a)') 'reading dataset: ', trim(dataset)

      call h5dopen_f(file_id, dataset, dset_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)')   'error: data set ',trim(dataset), ' could not be opened'
        write (*,'(a,i8)') 'error code ', ierror
        return
      endif

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

      if (maxndim.ne.max_dimen) then
         write (*,*) '  expecting ', max_dimen
         write (*,*) '  found     ', maxndim
         call h5_fatal('read_string1d','invalid number of dimensions')
      endif

!--- Note that the length of the character strings returned from the read
!--- are not the same as the lengths passed into the subroutine
!--- To fix, read data into a temporary character array then fill
!--- the actual character arrays

      if (h_dims(1).ne.nsize) then
         write (*,*) '  expecting ', nsize
         write (*,*) '  found     ', h_dims(1)
         call h5_fatal('read_string1d','invalid string size')
      endif

! ***** how can you determine klen???

      allocate (htemp(klen,nsize))

!--- determine data type

      call h5dget_type_f(dset_id, type_id, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dget_type_f ', ierror
        call h5_fatal('h5dget_type_f','ierror')
      endif

!  note that type_id does not match H5T_NATIVE_CHARACTER or H5T_STRING

!--- Read data from dataset - read temporary character array

      call h5dread_f(dset_id, type_id, htemp, h_dims, ierror)
      if (ierror.ne.0) then
        write (*,*) 'error: h5dread_f reading data set ',trim(dataset)
        call h5_fatal('read_string1d','error reading dataset')
      endif

! Close datatype and dataset

      call h5tclose_f(type_id, ierror)
      call h5dclose_f(dset_id, ierror)

! Transfer temporary character arrays to actual character arrays

      do j=1, nsize
        do i=1, klen
          if (ichar(htemp(i,j)).eq.0) then
            write (0,*) '*** found null character - fixing ****', i, j
            htemp(i,j)=' '   ! remove null  *** is this necessary?
          endif
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
!  Subroutine to write a scalar double
!
!  2017/08/08 - added, previously you had to write a scalar as a 1D array
!
!=======================================================================
      subroutine hwrite_double_scalar(file_id, dataset, xvar, update)
      use hdf5
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      real(8),          intent(in) :: xvar         ! scalar
      logical, optional,intent(in) :: update       ! set to true if update data instead of new data

!--- local

!!    integer,          intent(in) :: idims(10)    ! Dataset dimensions

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(hid_t) :: dtype_id      ! data type to be written
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

      integer  :: ierror              ! Error flag

      logical  :: ifupdate            ! local copy of optional argument

      if (present(update)) then
        ifupdate=update
      else
        ifupdate=.false.
      endif

!--- calculate rank from idims array

      dtype_id=H5T_NATIVE_DOUBLE

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset

      if (ifupdate) then     ! open existing dataset

        call h5dopen_f(file_id, dataset, dset_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5dopen_f','ierror')

      else

! Create a new dataspace (information about array)

!*** to create a true Scalar value (not a 1D array with one element),
!*** you could use "h5screate_f(H5S_SCALAR_F, dspace_id, ierror)   (not verified)

!!      call h5screate_f(dtype_id, dspace_id, ierror)
        call h5screate_f(H5S_SCALAR_F, dspace_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5screate_f','ierror')

! Create the dataset with default properties.

        call h5dcreate_f(file_id, dataset, dtype_id, dspace_id, dset_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5dcreate_f','ierror')

      endif

! Write dataset
!     (idimht is ignored for scalars)

      call h5dwrite_f(dset_id, dtype_id, xvar, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('h5dwrite_f','ierror')

! End access to the dataset and release resources used by it. (release dataset id)

      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5dclose_f','ierror')

! Terminate access to the data space.  (release dataspace id)

      if (.not. ifupdate) then
        call h5sclose_f(dspace_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5sclose_f','ierror')
      endif

      return
      end subroutine hwrite_double_scalar

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
!  8/3/2016 - updated with optional argument to allow subroutine to update
!             a dataset instead of creating a new one
!
!=======================================================================
      subroutine hwrite_double(file_id, dataset, idims, xvar, update)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: idims(10)    ! Dataset dimensions
      real(8),          intent(in) :: xvar(*)      ! multi-dimensional array, treated as 1D
      logical, optional,intent(in) :: update       ! set to true if update data instead of new data

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(hid_t) :: dtype_id      ! data type to be written

      integer  :: i
      integer  :: irank               ! Dataset rank
      integer  :: ierror              ! Error flag
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

      logical  :: ifupdate     ! local copy of optional argument

      if (present(update)) then
        ifupdate=update
      else
        ifupdate=.false.
      endif

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

      if (ifupdate) then     ! open existing dataset

        call h5dopen_f(file_id, dataset, dset_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5dopen_f','ierror')

      else

! Create a new dataspace (information about array)

!*** to create a true Scalar value (not a 1D array with one element),
!*** you could use "h5screate_f(H5S_SCALAR_F, dspace_id, ierror)   (not verified)

        call h5screate_simple_f(irank, idimht, dspace_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5screate_simple','ierror')

! Create the dataset with default properties.

        call h5dcreate_f(file_id, dataset, dtype_id, dspace_id, dset_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5dcreate_f','ierror')

      endif

! Write dataset

      call h5dwrite_f(dset_id, dtype_id, xvar, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('h5dwrite_f','ierror')

! End access to the dataset and release resources used by it. (release dataset id)

      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('h5dclose_f','ierror')

! Terminate access to the data space.  (release dataspace id)

      if (.not. ifupdate) then
        call h5sclose_f(dspace_id, ierror)
        if (ierror.ne.0) call h5_fatal('h5sclose_f','ierror')
      endif

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
      subroutine hwrite_integer_scalar(file_id, dataset, xvar)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: xvar         ! scalar

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(hid_t) :: dtype_id      ! data type to be written

      integer  :: ierror              ! Error flag
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

!--- calculate rank from idims array

      dtype_id=H5T_NATIVE_INTEGER

      write (*,'(2a,i3)') ' writing dataspace name: ', dataset

! Create the dataspace (information about scalar)

      call h5screate_f(H5S_SCALAR_F, dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5screate_f','ierror')

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
      end subroutine hwrite_integer_scalar

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
!=======================================================================
      subroutine hwrite_string1d(file_id, dataset, namex, idims)
      implicit none

      character(len=*), intent(in) :: dataset      ! Dataset name
      integer(hid_t),   intent(in) :: file_id      ! File identifier
      integer,          intent(in) :: idims        ! Dataset dimensions
      character(len=*), intent(in) :: namex(idims) ! 1D array

!--- local

      integer(hid_t) :: dset_id       ! Dataset identifier
      integer(hid_t) :: dspace_id     ! Dataspace identifier
      integer(size_t):: len_string    ! length of string
      integer(hid_t) :: type_id       ! attribute type to be written

      integer  :: ierror              ! Error flag
      integer  :: irank               ! dataset rank
      integer(hsize_t) :: idimht(10)  ! Dataset dimensions - using HDF size

! Create dataspace (information about array)

      irank=1               ! 1D array
      len_string=len(namex)
      idimht(1)=idims

!d    write (*,'(2a,i3)') ' writing dataspace name: ', dataset
!d    write (*,*)         ' debug: rank           : ', irank
!d    write (*,*)         ' debug: size           : ', idimht(1:irank)
!d    write (*,*)         ' debug: len            : ', len_string
      call h5screate_simple_f(irank, idimht, dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5screate_simple','ierror')

! Create data type id (from William)
      call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5tcopy_f','ierror')
      call h5tset_size_f(type_id, int(len_string, HSIZE_T), ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5tset_size_f','ierror')

! Create the dataset with default properties.

      call h5dcreate_f(file_id, dataset, type_id, dspace_id, dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dcreate_f','ierror')

! Create dataset attribute.

      call h5dwrite_f (dset_id, type_id, namex, idimht, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dwrite_f','ierror')

! Close the dataset

      call h5tclose_f(type_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5tclose_f','ierror')
      call h5dclose_f(dset_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5dclose_f','ierror')
      call h5sclose_f(dspace_id, ierror)
      if (ierror.ne.0) call h5_fatal('ierror: h5sclose_f','ierror')

      return
      end subroutine hwrite_string1d
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
