   program h_unit_test
!=======================================================================
!
!  Program to perform unit tests on Mod_htdftools
!
!  Program also serves as an example file to use Mod_hdftools routines
!
!  Copyright (c) 2014-2017 Core Physics, Inc.
!
!  Distributed under the MIT license.  
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2017/08/08 - add unit tests for reading/writing scalars
!
!=======================================================================

      use hdf5
      use mod_hdftools

      implicit none

      character(len=20) :: filename = 'dsetunit.h5'  ! File name
      character(len=12) :: dsetname                  ! Dataset name

      integer(hid_t)   :: file_id       ! File identifier

      integer     :: num1d
      integer     :: nx, ny, nz, na
      integer     :: i1, i2, i3, i4, i5
      integer     :: i, j
      real(8)     :: x1,xerr
      real(8)     :: tol=1.0d-8

      real(8)     :: power0d, atemp(1)
      real(8), allocatable :: power1d(:)
      real(8), allocatable :: power2d(:,:)
      real(8), allocatable :: power3d(:,:,:)
      real(8), allocatable :: power4d(:,:,:,:)
      real(8), allocatable :: power5d(:,:,:,:,:)

      integer     :: kerr
      integer     :: ierror              ! Error flag
      integer     :: nbad 
      integer     :: nfail               ! number of failures
      integer     :: ndim                ! number of dimensions in h5info
      integer     :: itype               ! data type in h5info
      integer     :: idims(10)           ! dataset dimensions
      integer     :: int0d
      integer     :: itemp(1)
      integer, allocatable :: int1d(:)
      integer, allocatable :: int2d(:,:)

      character(len=20) :: namex(3)
      character(len=22) :: stringin      ! make these strings different lengths
      character(len=18) :: stringout

!--- allocate and fill dummy arrays

      nfail=0

      num1d=44
      nx=4
      ny=6
      nz=2
      na=3

      power0d=3.14d0

      allocate (power1d(num1d))
      do i=1, num1d
        power1d(i)=dble(i*i)
      enddo

      allocate (power2d(nx,ny))
      do i2=1, ny
        do i1=1, nx
           power2d(i1,i2)=i1*100.0d0+i2+0.1d0
        enddo
      enddo

      allocate (power3d(nx,ny,nz))
      do i3=1, nz
        do i2=1, ny
          do i1=1, nx
            power3d(i1,i2,i3)=i3*1.0d3 + i1*100.0d0+i2+0.1d0
          enddo
        enddo
      enddo

      allocate (power4d(nx,ny,nz,na))
      do i4=1, na
        do i3=1, nz
          do i2=1, ny
            do i1=1, nx
              power4d(i1,i2,i3,i4)=i4*1.0d4 + i3*1.0d3 + i1*100.0d0 + i2+0.1d0
            enddo
          enddo
        enddo
      enddo

      allocate (power5d(nx,ny,nz,na,na))
      do i5=1, na
        do i4=1, na
          do i3=1, nz
            do i2=1, ny
              do i1=1, nx
                power5d(i1,i2,i3,i4,i5)=i5+1.0d5 + i4*1.0d4 + i3*1.0d3 + i1*100.0d0 + i2+0.1d0
              enddo
            enddo
          enddo
        enddo
      enddo

      int0d=3140

      allocate (int1d(num1d))
      do i=1, num1d
        int1d(i)=i*i
      enddo

      allocate (int2d(nx,ny))
      do j=1, ny
        do i=1, nx
           int2d(i,j)=i*100+j
        enddo
      enddo


!-------------------------------
!   write data to HDF file
!-------------------------------

!--- Initialize fortran interface.

      call h5open_f(ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

      write (*,'(2a,i10)') 'creating file: ', trim(filename)

!--- Create a new file using default properties.

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

!---------------

!--- call wrapper routine to write double scalar

      dsetname = 'power'
      call hwrite_double_scalar(file_id, dsetname, power0d)

!--- call wrapper routine to write double 0D array

      idims(:)=0     ! clear
      idims(1)=1
      atemp(1)=power0d   ! hack: send as single element array in order to pass type checking

      dsetname = 'power0d'
      call hwrite_double(file_id, dsetname, idims, atemp)

!--- call wrapper routine to write double 1D array

      idims(:)=0     ! clear
      idims(1)=num1d

      dsetname = 'power1d'
      call hwrite_double(file_id, dsetname, idims, power1d)

!--- call wrapper routine to write double 2D array

      idims(:)=0     ! clear
      idims(1)=nx
      idims(2)=ny

      dsetname = 'power2d'
      call hwrite_double(file_id, dsetname, idims, power2d)

!--- call wrapper routine to write double 3D array

      idims(:)=0     ! clear
      idims(1)=nx
      idims(2)=ny
      idims(3)=nz

      dsetname = 'power3d'
      call hwrite_double(file_id, dsetname, idims, power3d)

!--- call wrapper routine to write double 4D array

      idims(:)=0     ! clear
      idims(1)=nx
      idims(2)=ny
      idims(3)=nz
      idims(4)=na

      dsetname = 'power4d'
      call hwrite_double(file_id, dsetname, idims, power4d)

!--- call wrapper routine to write double 5D array

      idims(:)=0     ! clear
      idims(1)=nx
      idims(2)=ny
      idims(3)=nz
      idims(4)=na
      idims(5)=na

      dsetname = 'power5d'
      call hwrite_double(file_id, dsetname, idims, power5d)

!---------------

!--- call wrapper routine to write integer scalar

      dsetname = 'intscalar'
      call hwrite_integer_scalar(file_id, dsetname, int0d)

!--- call wrapper routine to write integer 0D array

      idims(:)=0     ! clear
      idims(1)=1
      itemp(1)=int0d   ! hack: send as single element array in order to pass type checking

      dsetname = 'int0d'
      call hwrite_integer(file_id, dsetname, idims, itemp)

!--- call wrapper routine to write integer 1D array

      idims(:)=0     ! clear
      idims(1)=num1d

      dsetname = 'int1d'
      call hwrite_integer(file_id, dsetname, idims, int1d)

!--- call wrapper routine to write integer 2D array

      idims(:)=0     ! clear
      idims(1)=nx
      idims(2)=ny

      dsetname = 'int2d'
      call hwrite_integer(file_id, dsetname, idims, int2d)

!---------------

!--- write string

      dsetname='bigdog'
      stringin='This dog can run'
      call hwrite_string(file_id, dsetname, stringin)

!--- write string array

      namex(1)='this bird'
      namex(2)='flew fast'
      namex(3)='over tree'

      idims(:)=0     ! clear
      idims(1)=3

! **** fix: this routine is not working, it only writes the first element of array ****

      dsetname='namex'
      call hwrite_string1d(file_id, dsetname, namex, idims(1))

!--- close file.

      call h5fclose_f(file_id, ierror)

      write (*,*) 'finished writing file'

!-----------------------------------
!   read data back from HDF file
!-----------------------------------

      write (*,*)

      write (*,'(2a)') 'reading h5 file: ', trim(filename)
      call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, ierror)
      if (ierror<0) then
        write(*,*) 'error: H5 input file ',trim(filename),' could not be opened'
        stop
      endif

!--- test h5info

      dsetname='junk'
      call h5info(file_id, dsetname, itype, ndim, idims)  ! test return code for not found
      if (itype.eq.-99) then
        write (*,220) 'h5info return code 1 ', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'h5finfo did not return proper return code'
        write (*,220) 'h5info return code 1   ', 'FAIL'
      endif

      dsetname='power4d'
      call h5info(file_id, dsetname, itype, ndim, idims)

      if (itype.eq.0) then
        write (*,220) 'h5info return code 2   ', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'h5finfo did not return correct itype'
        write (*,220) 'h5info return code 2   ', 'FAIL'
      endif

      if (ndim.eq.4) then
        write (*,220) 'h5info number of dimensions', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'h5finfo returned invalid number of dimensions'
        write (*,220) 'h5info number of dimensions', 'FAIL'
      endif

      call checkint(1, nx,   idims(1), nfail)
      call checkint(2, ny,   idims(2), nfail)
      call checkint(3, nz,   idims(3), nfail)
      call checkint(4, na,   idims(4), nfail)
      call checkint(5, 0,    idims(5), nfail)

!--- read scalar double

      dsetname = 'power'
      power0d=0.0d0        ! clear before read

      call hdf5_read_double (file_id, dsetname, power0d)
      xerr=abs(power0d-3.14d0)
      if (xerr.lt.tol) then
        write (*,220) 'read_double', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'value returned ', power0d
        write (*,*) 'expecting      ', 3.14d0
        write (*,*) 'error = ', xerr
        write (*,220) 'read_double', 'FAIL'
      endif

!--- read 0D array of double

      dsetname = 'power0d'
      power0d=0.0d0        ! clear before read

      call hdf5_read_double (file_id, dsetname, power0d)
      xerr=abs(power0d-3.14d0)
      if (xerr.lt.tol) then
        write (*,220) 'read_double', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'value returned ', power0d
        write (*,*) 'expecting      ', 3.14d0
        write (*,*) 'error = ', xerr
        write (*,220) 'read_double', 'FAIL'
      endif

!--- read 1D array of double

      dsetname = 'power1d'
      power1d=0.0d0        ! clear before read

      call hdf5_read_double (file_id, dsetname, num1d, i, power1d)
      if (i.ne.num1d) then
         write (*,*) 'read 1D array size ', i
         write (*,*) 'expecting          ', num1d
         write (*,220) 'error in returned 1D size', 'FAIL'
         nfail=nfail+1
      endif
      nbad=0
      do i=1, num1d
        xerr=abs(power1d(i)-dble(i*i))
        if (xerr.gt.tol) nbad=nbad+1  
      enddo
      call printstatus('read_double1d', nbad, nfail)

!--- read 2D array of double

      dsetname = 'power2d'
      power2d=0.0d0

      call hdf5_read_double (file_id, dsetname, nx, ny, power2d)
      nbad=0
      do j=1, ny
        do i=1, nx
          xerr=abs(power2d(i,j)-(i*100.0d0+j+0.1d0))
          if (xerr.gt.tol) nbad=nbad+1
        enddo
      enddo
      call printstatus('read_double2d', nbad, nfail)

!--- read 3D array of double

      dsetname = 'power3d'
      power3d=0.0d0

      call hdf5_read_double (file_id, dsetname, nx, ny, nz, power3d)
      nbad=0
      do i3=1, nz
        do i2=1, ny
          do i1=1, nx
            xerr=abs(power3d(i1,i2,i3)-(i3*1000.0d0+i1*100.0d0+i2+0.1d0))
            if (xerr.gt.tol) nbad=nbad+1
           enddo
        enddo
      enddo
      call printstatus('read_double3d', nbad, nfail)

!--- read 4D array of double

      dsetname = 'power4d'
      power4d=0.0d0
      
      call hdf5_read_double (file_id, dsetname, nx, ny, nz, na, power4d)
      nbad=0
      do i4=1, na
        do i3=1, nz
          do i2=1, ny
            do i1=1, nx
              x1=i4*10000.0d0+i3*1000.0d0+i1*100.0d0+i2+0.1d0
              xerr=abs(power4d(i1,i2,i3,i4)-x1)
              if (xerr.gt.tol) nbad=nbad+1
            enddo
          enddo
        enddo
      enddo
      call printstatus('read_double4d', nbad, nfail)

!--- read 5D array of double

      dsetname = 'power5d'
      power5d=0.0d0

      call hdf5_read_double (file_id, dsetname, nx, ny, nz, na, na, power5d)
      nbad=0
      do i5=1, na
        do i4=1, na
          do i3=1, nz
            do i2=1, ny
              do i1=1, nx
                x1=i5+1.0d5+i4*1.0d4+i3*1000.0d0+i1*100.0d0+i2+0.1d0
                xerr=abs(power5d(i1,i2,i3,i4,i5)-x1)
                if (xerr.gt.tol) nbad=nbad+1
              enddo
            enddo
          enddo
        enddo
      enddo
      call printstatus('read_double5d', nbad, nfail)

!--- read scalar integer

      dsetname = 'intscalar'
      int0d=0        ! clear before read

      call hdf5_read_integer (file_id, dsetname, int0d)
      kerr=int0d-3140
      if (kerr.eq.0) then
        write (*,220) 'read_integer_scalar', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'value returned ', int0d
        write (*,*) 'expecting      ', 3140
        write (*,*) 'error = ', kerr
        write (*,220) 'read_integer', 'FAIL'
      endif


!--- read 0D array of integer

      dsetname = 'int0d'
      int0d=0        ! clear before read

      call hdf5_read_integer (file_id, dsetname, int0d)
      kerr=int0d-3140
      if (kerr.eq.0) then
        write (*,220) 'read_integer', 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'value returned ', int0d
        write (*,*) 'expecting      ', 3140
        write (*,*) 'error = ', kerr
        write (*,220) 'read_integer', 'FAIL'
      endif

!--- read 1D array of integer (max size used)

      dsetname = 'int1d'
      int1d=0        ! clear before read

      call hdf5_read_integer (file_id, dsetname, num1d, i, int1d)
      if (i.ne.num1d) then
         write (*,*) 'read 1D array size ', i
         write (*,*) 'expecting          ', num1d
         write (*,220) 'error in returned 1D size', 'FAIL'
         nfail=nfail+1
      endif
      nbad=0
      do i=1, num1d
        kerr=int1d(i)-i*i
        if (kerr.ne.0) nbad=nbad+1
      enddo
      call printstatus('read_integer1d', nbad, nfail)

!--- read 1D array of integer (size known)

      dsetname = 'int1d'
      int1d=0        ! clear before read

      call hdf5_read_integer (file_id, dsetname, num1d, int1d)
      nbad=0
      do i=1, num1d
        kerr=int1d(i)-i*i
        if (kerr.ne.0) nbad=nbad+1
      enddo
      call printstatus('read_integer1d', nbad, nfail)

!--- read 2D array of integer

      dsetname = 'int2d'
      int2d=0

      call hdf5_read_integer (file_id, dsetname, nx, ny, int2d)
      nbad=0
      do j=1, ny
        do i=1, nx
          kerr=int2d(i,j)-(i*100+j)
          if (kerr.ne.0) nbad=nbad+1
        enddo
      enddo
      call printstatus('read_integer2d', nbad, nfail)

!--- read string

      dsetname = 'bigdog'
      stringout='X'

      call read_string(file_id, dsetname, stringout)

      if (stringout.eq.stringin) then
        write (*,220) 'read_string', 'PASS'
      else
        nfail=nfail+1
        write (*,220) 'read_string', 'FAIL'
      endif

!--- read 1D array of strings

      namex(1)='X'
      namex(2)='X'
      namex(3)='X'
      i=3   ! pass in number to read

      dsetname='namex'
      call read_string1d(file_id, dsetname, namex, i)

      nbad=0
      if (namex(1).ne.'this bird') nbad=nbad+1
      if (namex(2).ne.'flew fast') nbad=nbad+1
      if (namex(3).ne.'over tree') nbad=nbad+1
      if (nbad.eq.0) then
        write (*,220) 'read_string1d', 'PASS'
      else
        write (*,*) 'number of errors ', nbad
        nfail=nfail+1
        write (*,220) 'read_string1d', 'FAIL'
      endif

!--- close file

      call h5fclose_f(file_id, ierror)

!-----------------------------------
!   Test ability to overwrite data
!-----------------------------------

      write (*,*)
      write (*,*)  'Testing ability to overwrite existing data'

      call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, ierror)
      if (ierror<0) then
        write(*,*) 'error: H5 input file ',trim(filename),' could not be opened'
        stop
      endif

!--- overwrite double 2D array

      do j=1, ny
        do i=1, nx
           power2d(i,j)=i*200.0d0+j
        enddo
      enddo

      idims(:)=0     ! clear
      idims(1)=nx
      idims(2)=ny

      dsetname = 'power2d'
      call hwrite_double(file_id, dsetname, idims, power2d, update=.true.)

!--- read back

      dsetname = 'power2d'
      power2d=0.0d0

      call hdf5_read_double (file_id, dsetname, nx, ny, power2d)
      nbad=0
      do j=1, ny
        do i=1, nx
          xerr=abs(power2d(i,j)-(i*200.0d0+j))
          if (xerr.gt.tol) nbad=nbad+1
        enddo
      enddo
      call printstatus('read_double2d', nbad, nfail)

!--- close file

      call h5fclose_f(file_id, ierror)


!-------------------------------
!   finished
!-------------------------------

  220 format(1x,a,t60,a)

!--- Close fortran interface.

      call h5close_f(ierror)

      deallocate (power1d)
      deallocate (power2d)

!--- if no failures, delete temporary file

      if (nfail.eq.0) then
        write (*,*)
        write (*,*) 'deleting temporary h5 file since everything passed'
        open (20,file='dsetunit.h5', status='old')
        close (20,status='delete')
      endif

!--- print final results

      write (*,*)
      if (nfail.eq.0) then
        write (*,220) 'Overall', 'PASS'
      else
        write (*,220) 'Overall', 'FAIL'
      endif


      end
!=======================================================================
!  Small subroutine to check two integer values
!=======================================================================
      subroutine checkint(id, i1,i2, nfail)
      implicit none
      integer, intent(in) :: id      ! dimension number (edit)
      integer, intent(in) :: i1, i2  ! expected, actual
      integer             :: nfail   ! return code

      if (i1.eq.i2) then
        write (*,220) id, 'PASS'
      else
        nfail=nfail+1
        write (*,*) 'h5finfo returned invalid dimension size ', id
        write (*,*) ' expecting ', i1
        write (*,*) ' found     ', i2
        write (*,220) id, 'FAIL'
      endif
  220 format(' h5info dimension ', i0, t60,a)

      return
      end subroutine checkint

!=======================================================================
!  Small subroutine to print error messages
!=======================================================================
      subroutine printstatus(label, nbad, nfail)
      implicit none
      character(len=*), intent(in) :: label
      integer,          intent(in) :: nbad    ! current test
      integer                      :: nfail   ! running total

      if (nbad.eq.0) then
        write (*,220) trim(label), 'PASS'
      else
        write (*,*) 'number of errors ', nbad
        nfail=nfail+1
        write (*,220) trim(label), 'FAIL'
      endif
  220 format(1x,a,t60,a)

      return
      end subroutine printstatus

