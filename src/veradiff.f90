    program veradiff
!=======================================================================
!
!  Program to compare twoH VERA HDF output file and print summary
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
!  This program is intended to be used in regressing testing
!
!  Right now, we are just performing a simple diff on eigenvalues with a tolerance of 10 pcm.
!  This will be modifed in the future to check other parameters as well.
!
!  2013/09/06 - original file
!  2014/01/27 - compare boron if dataset is present on file
!  2014/04/26 - add comparison for multiple statepoints - see #3164
!  2014/11/01 - added pin power comparisons             - see #3430
!
!
!--------------------------------------------------------------------------------
      use hdf5
      implicit none

      character(len=250) :: fname1, fname2   ! HDF file names
      integer            :: iargs            ! number of command line arguments

      integer(hid_t)     :: file_id1, file_id2  ! HDF file id's
      integer            :: ierror
      integer            :: nfail            ! number of diff failures
      integer            :: nstate           ! statepoint number

      logical            :: ifxst

      character(len=10) :: group_name

      ierror=0
      fname1=' '
      fname2=' '

      nfail=0    ! number of test failures

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.2) then
        write (*,*) 'usage:  veradiff.exe [file1] [file2]'
        nfail=nfail+1
        goto 900
      endif

      call get_command_argument(1,fname1)
      call get_command_argument(2,fname2)

      write (*,'(/,1x,2a)') 'reading h5 file: ', trim(fname1)
      inquire(file=fname1, exist=ifxst)
      if (.not.ifxst) then
        write (*,*) 'ERROR: input file ',trim(fname1),' does not exist'
        nfail=nfail+1
        goto 900
      endif

      write (*,'(1x,2a)') 'reading h5 file: ', trim(fname2)
      inquire(file=fname2, exist=ifxst)
      if (.not.ifxst) then
        write (*,*) 'ERROR: input file ',trim(fname2),' does not exist'
        nfail=nfail+1
        goto 900
      endif

!--------------------------------------------------------------------------------
! Initialize HDF and open files
!--------------------------------------------------------------------------------

!--- initialize fortran interface

      call h5open_f(ierror)      ! NOTE: THIS IS REQUIRED

!--- open file1

      call h5fopen_f (fname1, H5F_ACC_RDWR_F, file_id1, ierror)
      if (ierror<0) then
        write (*,'(3a)') 'ERROR: H5 input file ',trim(fname1),' could not be opened'
        nfail=nfail+1
        goto 900
      endif

!--- open file2

      call h5fopen_f (fname2, H5F_ACC_RDWR_F, file_id2, ierror)
      if (ierror<0) then
        write (*,'(3a)') 'ERROR: H5 input file ',trim(fname2),' could not be opened'
        nfail=nfail+1
        goto 900
      endif

!--- check if "STATE_0001" exists on file.  If so, this is a new file

      group_name='STATE_0001'
      call h5lexists_f(file_id1, group_name, ifxst, ierror)
      if (ifxst) then
        nstate=1
      else
        nstate=0    ! old file without STATE groups
      endif

!--- loop over all statepoints and do comparisons

      do
        call readstate(file_id1, file_id2, nstate, nfail)
        if (nstate.eq.0) exit   ! old files only have one statepoint

        nstate=nstate+1
        write (group_name(7:10),'(i4.4)') nstate
        call h5lexists_f(file_id1, group_name, ifxst, ierror)
        if (.not.ifxst) exit

      enddo

!--------------------------------------------------------------------------------
! Close File
!--------------------------------------------------------------------------------
      call h5fclose_f(file_id1, ierror)
      call h5fclose_f(file_id2, ierror)

!--------------------------------------------------------------------------------
! Print Results
!--------------------------------------------------------------------------------

  900 continue

      if (nfail.gt.0) then
        write (*,'(/,a)') 'Overall FAIL'
        call exit(nfail)    ! return error code
      else
        write (*,'(/,a)') 'Overall PASS'
      endif

      end program

!=======================================================================
!
!  Subroutine to read a single statepoint from HDF file
!
!=======================================================================
      subroutine readstate(file_id1, file_id2, nstate, nfail)
      use hdf5
      use mod_hdftools, only : hdf5_read_double
      implicit none

      integer(hid_t), intent(in) :: file_id1, file_id2
      integer,        intent(in) :: nstate
      integer                    :: nfail

!--- local variables

      integer            :: ierror

      logical            :: ifmissing
      logical            :: ifxst

      real(8)            :: boron1, boron2
      real(8)            :: bordiff
      real(8)            :: bortol=1.0d0   ! boron tolerance (ppm)

      real(8)            :: xkeff1, xkeff2
      real(8)            :: xkdiff
      real(8)            :: xktol=10.0d0   ! eigenvalue tolerance (pcm)

      real(8)            :: rms

      character(len=80)  :: dataset
      character(len=12)  :: statename

      boron1=-9999.0d0    ! large negative number if dataset not present
      boron2=-9999.0d0    ! large negative number if dataset not present

      xkeff1=0.0d0
      xkeff2=0.0d0

      ierror=0
      ifmissing=.false.

      write (*,*)

      if (nstate.eq.0) then
        statename=' '
      else
        statename='STATE_0000/'
        write (statename(7:10),'(i4.4)') nstate
      endif

!--------------------
!  eigenvalue
!--------------------

      dataset=trim(statename)//'keff'

      call h5lexists_f(file_id1, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id1, dataset, xkeff1)
      else
        nfail=nfail+1
        ifmissing=.true.
        write (*,*) 'FAIL - keff dataset does not exist on file 1'
      endif

      call h5lexists_f(file_id2, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id2, dataset, xkeff2)
      else
        nfail=nfail+1
        ifmissing=.true.
        write (*,*) 'FAIL - keff dataset does not exist on file 2'
      endif

      if (ifmissing) return    ! return if eigenvalue missing

      xkdiff=(xkeff2-xkeff1)*1.0d5

      write (*,'(/,1x,a)') 'Eigenvalue comparisons'
      write (*,'(2x,a, f12.7)')   'keff1  =', xkeff1
      write (*,'(2x,a, f12.7)')   'keff2  =', xkeff2
      write (*,'(2x,a, f12.2,a)') 'diff   =', xkdiff,' pcm'
      write (*,'(2x,a, f12.2,a)') 'tol    =', xktol ,' pcm'

      if (abs(xkdiff).lt.xktol) then
        write (*,*) 'PASS - eigenvalue difference is less than tolerance'
      else
        nfail=nfail+1
        write (*,*) 'FAIL - eigenvalue difference exceeds tolerance'
      endif

!--------------------
!  boron (optional)
!--------------------

      write (*,*)

      dataset=trim(statename)//'boron'

      call h5lexists_f(file_id1, dataset, ifxst, ierror)
      if (ifxst) call hdf5_read_double(file_id1, dataset, boron1)
      call h5lexists_f(file_id2, dataset, ifxst, ierror)
      if (ifxst) call hdf5_read_double(file_id2, dataset, boron2)

!  compare boron

      if (boron1.gt.-1000.0d0 .and. boron2.gt.-1000.0d0) then
        write (*,'(/,1x,a)') 'Boron comparisons'
        bordiff=(boron2-boron1)
        write (*,'(2x,a, f12.2)')   'boron1 =', boron1
        write (*,'(2x,a, f12.2)')   'boron2 =', boron2
        write (*,'(2x,a, f12.2,a)') 'diff   =', bordiff,' ppm'
        write (*,'(2x,a, f12.2,a)') 'tol    =', bortol ,' ppm'
        if (abs(bordiff).lt.bortol) then
          write (*,*) 'PASS - boron difference is less than tolerance'
        else
          nfail=nfail+1
          write (*,*) 'FAIL - boron difference exceeds tolerance'
        endif
      else
        write (*,'(/,1x,a)') 'No boron comparison performed - at least one file is missing boron'
      endif

!----------------
!  pin powers
!----------------

      dataset=trim(statename)//'pin_powers'

      call h5lexists_f(file_id1, dataset, ifxst, ierror)
      if (.not.ifxst) then
        nfail=nfail+1
        ifmissing=.true.
        write (*,*) 'FAIL - pin power dataset does not exist on file 1'
      endif

      call h5lexists_f(file_id2, dataset, ifxst, ierror)
      if (.not.ifxst) then
        nfail=nfail+1
        ifmissing=.true.
        write (*,*) 'FAIL - pin power dataset does not exist on file 2'
      endif

      if (.not.ifmissing) then
        call pin_compare(file_id1, file_id2, dataset, rms)
        write (*,*) 'No PIN tolerance used'
      endif

! ****** NOTE THAT NO TOLERANCE IS SET, PIN POWER DIFFERENCES *******
! ******        IS JUST INFORMATIONAL AT THIS POINT           *******

!--- return

      return

      end subroutine readstate
!=======================================================================
!
!  Subroutine to read two pin_power datasets and calculate RMS
!
!=======================================================================
      subroutine pin_compare(file_id1, file_id2, dataset, rms)
      use hdf5
      use mod_hdftools, only : h5info, hdf5_read_double
      implicit none
      integer(hid_t),   intent(in)  :: file_id1, file_id2
      character(len=*), intent(in)  :: dataset
      real(8),          intent(out) :: rms

      integer :: ip, jp, k, na, np
      integer :: nassm, kd, npin
      integer :: ndim         ! number of array dimensions
      integer :: idim(10)     ! array dimension sizes
      integer :: itype        !

      logical :: ifbad        ! flag for bad data

      real(8)              :: p1, dd, dmax
      real(8), allocatable :: power1(:,:,:,:)
      real(8), allocatable :: power2(:,:,:,:)

!  find dimensions so we can allocate

      rms=999999.0d0   ! set high in case reads fail
      ifbad=.false.

      call h5info(file_id1, dataset, itype, ndim, idim)

      nassm=idim(1)    ! order (nassm, kd, npin, npin)
      kd   =idim(2)
      npin =idim(4)

!d    write (*,*) 'debug pin: ndim  = ', ndim
!d    write (*,*) 'debug pin: nassm = ', nassm
!d    write (*,*) 'debug pin: kd    = ', kd
!d    write (*,*) 'debug pin: npin  = ', npin

      if (ndim .ne.4) then
        ifbad=.true.
        write (*,*) 'invalid dimensions in pin data'
      endif
      if (nassm.eq.0) then
        ifbad=.true.
        write (*,*) 'invalid number of assemblies in pin data'
      endif
      if (npin .eq.0) then
        ifbad=.true.
        write (*,*) 'invalid number of pins in pin data'
      endif
      if (idim(3).ne.idim(4)) then
        ifbad=.true.
        write (*,*) 'invalid npin in pin data'
      endif
      if (ifbad) return   ! return with big RMS

      allocate (power1(nassm, kd, npin, npin))
      allocate (power2(nassm, kd, npin, npin))

      call hdf5_read_double(file_id1, dataset, idim(1), idim(2), idim(3), idim(4), power1)
      call hdf5_read_double(file_id2, dataset, idim(1), idim(2), idim(3), idim(4), power2)

      rms=0.0d0
      dmax=0.0d0
      np=0
      do ip=1, npin
        do jp=1, npin
          do k=1, kd
            do na=1, nassm
              p1=power1(na,k,jp,ip)
              if (p1.gt.0.0d0) then
                np=np+1
                dd=power2(na,k,jp,ip)-p1
                rms=rms+dd*dd
                dmax=max(dmax,dd)
              endif
            enddo
          enddo
        enddo
      enddo

      if (np.gt.0) rms=sqrt(rms/dble(np))
 
      deallocate (power2)
      deallocate (power1)

      write (*,'(/1x,a)') 'Pin Power comparisons'
      write (*,*)  ' num values=', np
      write (*,30) 'max', dmax*100.0d0
      write (*,30) 'rms', rms*100.0d0
  30  format (2x,a,' difference =', f10.4,' %')

      return
      end subroutine pin_compare
!=======================================================================
