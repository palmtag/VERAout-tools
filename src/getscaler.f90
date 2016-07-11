   program getscaler
!=======================================================================
!
!  Program to read a single (double precision) scalar value from HDF file
!
!  *** needs to be modifed to read any data type ***
!
!  Copyright (c) 2016 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2016/06/21 - initial draft
!
!  ****** currently only works for type double *****
!  ****** currently only works for type double *****
!  ****** currently only works for type double *****
!
!  Usage:
!    getscaler [dataset] [file.h5] {reference} {reference tolerance}
!
!  Example:
!    getscaler /STATE_0001/tinlet tin3.h5 291 0.1
!
!-----------------------------------------------------------------------
      use  hdf5
      use  mod_hdftools
      implicit none

      character(len=80)  :: inputfile
      character(len=80)  :: reference       ! reference value
      character(len=80)  :: reftol          ! reference tolerance
      character(len=80)  :: dataset         ! HDF dataset name
      character(len=4)   :: result          ! pass or fail

      integer            :: iargs           ! number of command line arguments
      integer            :: ierror

      real(8)            :: xtmp
      real(8)            :: xdiff
      real(8)            :: reftmp    ! reference
      real(8)            :: tol       ! tolerance

      logical            :: ifxst

      integer(hid_t)     :: file_id

!  initialize

      inputfile=' '
      dataset  =' '
      reference=' '
      reftol   =' '
      result   =' '

!--- Read in arguments from command line

      iargs = command_argument_count()
      if (iargs.lt.2) then
        write (*,*) ' usage:'
        write (*,*) 'getscaler [dataset] [file.h5] {reference} {reference tolerance}'
        stop
      endif

      call get_command_argument(1,dataset)
      call get_command_argument(2,inputfile)

      if (iargs.ge.3) call get_command_argument(3,reference)
      if (iargs.ge.4) call get_command_argument(4,reftol)

      if (inputfile.eq.' ') then
        stop 'no input file specified on command line'
      endif

!--- initialize fortran interface

      call h5open_f(ierror)
      if (ierror.lt.0) then
        write (*,'(a)') 'error opening HDF file'
        stop
      endif

!--- open HDF file

      write (*,'(2a)') ' reading h5 file: ', trim(inputfile)
      inquire(file=inputfile, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(inputfile),' does not exist'
        stop
      endif
      call h5fopen_f (inputfile, H5F_ACC_RDONLY_F, file_id, ierror)   ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop
      endif

!--- check if dataset exists

      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (.not.ifxst) then
        write (0,*) 'dataset does not exist on file'
        stop
      endif

!*********** NEED TO CHECK DATASET TYPE ******************
!*********** NEED TO CHECK DATASET TYPE ******************
!*********** NEED TO CHECK DATASET TYPE ******************
!*********** NEED TO CHECK DATASET TYPE ******************

!--- read dataset

      call hdf5_read_double(file_id, dataset, xtmp)

!!    call hdf5_read_integer(file_id, dataset, itmp)
!!    call read_string(file_id, dataset, title)

!--- close file

! 800 continue

      call h5fclose_f(file_id, ierror)

!--- compare to reference

      if (reference.ne.' ') then
         read (reference,*) reftmp
         write (*,*) 'reference = ', reftmp
         xdiff=xtmp-reftmp

         if (reftol.eq.' ') then
            tol=0.01d0
         else
            read (reftol,*) tol
         endif
         write (*,*) 'reftol = ', tol

         if (abs(xdiff).lt.tol) then
            result='PASS'
         else
            result='FAIL'
         endif

      endif

!--- print summary

      write (*,*)
      write (*,*) xtmp
      if (result.ne.' ') write (*,*) trim(result)

!--- finished

      end program

