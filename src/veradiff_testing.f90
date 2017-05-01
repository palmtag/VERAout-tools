   program veradiff
   use hdf5
   implicit none
!=======================================================================
!
!  Program to compare two VERA HDF output file and print summary
!
!  **** This version of the diff program is intended to be used  ****
!  **** with automated testing where you need a pass/fail result ****
!
!  Copyright (c) 2014-2016 Core Physics, Inc.
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
!  2014/12/05 - added pin power tolerance               - see #3430
!  2016/02/28 - add statepoint summary
!
!--------------------------------------------------------------------------------

      character(len=250) :: fname1, fname2   ! HDF file names
      integer            :: iargs            ! number of command line arguments

      integer(hid_t)     :: file_id1, file_id2  ! HDF file id's
      integer            :: ierror
      integer            :: nfail            ! number of diff failures
      integer            :: n
      integer            :: nstate           ! statepoint number

      logical            :: ifxst

      real(8), allocatable :: summary(:,:)

      ierror=0
      fname1=' '
      fname2=' '

      nstate=0   ! number of statepoints on file
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
      endif

      write (*,'(1x,2a)') 'reading h5 file: ', trim(fname2)
      inquire(file=fname2, exist=ifxst)
      if (.not.ifxst) then
        write (*,*) 'ERROR: input file ',trim(fname2),' does not exist'
        nfail=nfail+1
      endif

      if (nfail.gt.0) goto 900

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

!--- count number of statepoints

      call countstates(file_id1, nstate)
      if (nstate.eq.0) nfail=nfail+1

      call countstates(file_id2, n)
      if (n.eq.0) nfail=nfail+1

      if (nstate.ne.n) then
        write (*,*) 'Number of statepoints on file 1 = ', nstate
        write (*,*) 'Number of statepoints on file 2 = ', n
        write (*,'(a)') 'ERROR: Number of statepoints do not match'
        nstate=0      ! protect summary messages
        nfail=nfail+1
      endif

      if (nfail.gt.0) goto 900

      allocate (summary(7,nstate))
      summary=0.0d0

!--- loop over all statepoints and do comparisons

      do n=1, nstate
        call readstate(file_id1, file_id2, n, nfail, summary(1,n))
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

      if (nstate.gt.1) then
        write (*,122)
        do n=1, nstate
          write (*,55) n, summary(:,n)
        enddo
      endif
  122 format ( &
          /,'  ==================================================================', &
          /,'                     Statepoint Summary (file2-file1)', &
          /,'  ==================================================================', &
          /,'     N     exp     EFPD    keff1     keff2      pcm   max%    rms%')
   55 format (2x,i4, f9.4, f8.2, 2f10.6, f8.1, 2f8.4)

      if (allocated(summary)) deallocate (summary)

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
      subroutine readstate(file_id1, file_id2, nstate, nfail, summary)
      use hdf5
      use mod_hdftools, only : hdf5_read_double
      implicit none

      integer(hid_t), intent(in) :: file_id1, file_id2
      integer,        intent(in) :: nstate
      integer                    :: nfail
      real(8)                    :: summary(7)    ! summary output

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

      real(8)            :: xexpo          ! exposure
      real(8)            :: xefpd          ! exposure in EFPD
      real(8)            :: xmax           ! max pin difference (%)
      real(8)            :: xrms           ! RMS pin difference (%)
      real(8)            :: pintol=0.5d0   ! pin power tolerance (%)

      character(len=80)  :: dataset
      character(len=12)  :: statename

      boron1=-9999.0d0    ! large negative number if dataset not present
      boron2=-9999.0d0    ! large negative number if dataset not present

      xkeff1=0.0d0
      xkeff2=0.0d0

      ierror=0
      ifmissing=.false.

      write (*,*)

      statename='STATE_0000/'
      write (statename(7:10),'(i4.4)') nstate

      write (*,120) nstate
  120 format (' =======================================', &
            /,'       Statepoint ', i0, &
            /,' =======================================')

!--------------------
!  exposure
!--------------------

      dataset=trim(statename)//'exposure'
      call hdf5_read_double(file_id1, dataset, xexpo)

      dataset=trim(statename)//'exposure_efpd'
      call hdf5_read_double(file_id1, dataset, xefpd)

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
        call pin_compare(file_id1, file_id2, dataset, xmax, xrms)
        write (*,320) pintol
        if (abs(xmax).lt.pintol) then
          write (*,*) 'PASS - pin power difference is less than tolerance'
        else
          nfail=nfail+1
          write (*,*) 'FAIL - pin power difference exceeds tolerance'
        endif
      endif
  320 format (3x,'tolerance = ', f5.2,' %')

!--- print summary line

      summary(1)=xexpo
      summary(2)=xefpd
      summary(3)=xkeff1
      summary(4)=xkeff2
      summary(5)=xkdiff
      summary(6)=xmax
      summary(7)=xrms

      write (*,55) xkeff1, xkeff2, xkdiff, xmax, xrms
   55 format ('Summary  ',2f10.6, f8.1, 2f8.4)

!--- return

      return

      end subroutine readstate
!=======================================================================
!
!  Subroutine to read two pin_power datasets and calculate RMS
!
!=======================================================================
      subroutine pin_compare(file_id1, file_id2, dataset, xmax, xrms)
      use hdf5
      use mod_hdftools, only : h5info, hdf5_read_double
      implicit none
      integer(hid_t),   intent(in)  :: file_id1, file_id2
      character(len=*), intent(in)  :: dataset
      real(8),          intent(out) :: xmax   ! return max difference (%)
      real(8),          intent(out) :: xrms   ! return RMS difference (%)

      integer :: ip, jp, k, na, np
      integer :: nassm, kd, npin
      integer :: ndim         ! number of array dimensions
      integer :: idim(10)     ! array dimension sizes
      integer :: itype        !
      integer :: kmax1(4)     ! location of max difference
      integer :: lmax1(4)     ! location of max power file 1
      integer :: lmax2(4)     ! location of max power file 2

      logical :: ifbad        ! flag for bad data

      real(8)              :: p1, p2, dd, dmax
      real(8)              :: pmax1, pave1
      real(8)              :: pmax2, pave2
      real(8)              :: rms
      real(8), allocatable :: power1(:,:,:,:)
      real(8), allocatable :: power2(:,:,:,:)

!  find dimensions so we can allocate

      xmax=999999.0d0   ! set high in case reads fail
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
      if (ifbad) return   ! return with big xmax

      allocate (power1(nassm, kd, npin, npin))
      allocate (power2(nassm, kd, npin, npin))

      call hdf5_read_double(file_id1, dataset, idim(1), idim(2), idim(3), idim(4), power1)
      call hdf5_read_double(file_id2, dataset, idim(1), idim(2), idim(3), idim(4), power2)

      rms=0.0d0
      dmax=0.0d0
      np=0
      kmax1(:)=0    ! location of max difference
      lmax1(:)=0    ! location of max power file 1
      lmax2(:)=0    ! location of max power file 2
      pmax1=0.0d0   ! max power of distribution 1
      pmax2=0.0d0   ! max power of distribution 2
      pave1=0.0d0   ! ave power of distribution 1
      pave2=0.0d0   ! ave power of distribution 2
      do ip=1, npin
        do jp=1, npin
          do k=1, kd
            do na=1, nassm
              p1=power1(na,k,jp,ip)
              p2=power2(na,k,jp,ip)
              if (p1.gt.0.0d0) then
                np=np+1
                pave1=pave1+p1
                pave2=pave2+p2
                if (p1.gt.pmax1) then
                  pmax1=p1
                  lmax1(1)=ip
                  lmax1(2)=jp
                  lmax1(3)=k
                  lmax1(4)=na
                endif
                if (p2.gt.pmax2) then
                  pmax2=p2
                  lmax2(1)=ip
                  lmax2(2)=jp
                  lmax2(3)=k
                  lmax2(4)=na
                endif
                dd=abs(p2-p1)
                rms=rms+dd*dd
                if (dd.gt.dmax) then
                  dmax=dd
                  kmax1(1)=ip
                  kmax1(2)=jp
                  kmax1(3)=k
                  kmax1(4)=na
                endif
              endif
            enddo
          enddo
        enddo
      enddo

      if (np.gt.0) then
        rms=sqrt(rms/dble(np))
        pave1=pave1/dble(np)
        pave2=pave2/dble(np)
      endif

      xmax=dmax*100.0d0
      xrms=rms*100.0d0
 
      deallocate (power2)
      deallocate (power1)

      write (*,'(/1x,a)') 'Pin Power comparisons'
      write (*,120) np
      write (*,126) 1, pmax1, lmax1(:)
      write (*,126) 2, pmax2, lmax2(:)
      write (*,125) 1, pave1
      write (*,125) 2, pave2
      write (*,131) dmax*100.0d0, kmax1(:)
      write (*,130) rms*100.0d0
      write (*,138) (pmax2-pmax1)*100.0d0
 120  format (3x,'number values',2x,i10)
 125  format (3x,'ave power',i1,5x, f10.4,'   (does not include volume weighting)')
 126  format (3x,'max power',i1,5x, f10.4,'   at (ip,jp,k,na) ', 4i4)
 130  format (3x,'rms difference ', f10.4,' %')
 131  format (3x,'max difference ', f10.4,' % at (ip,jp,k,na) ', 4i4)
 138  format (3x,'difference in max pin ', f10.4,' %')


      return
      end subroutine pin_compare
!=======================================================================
!
!  Subroutine to count number of statepoints on a file
!  This is needed to perform allocations
!
!=======================================================================
      subroutine countstates(file_id, nstate)
      use hdf5
      implicit none
      integer(hid_t), intent(in)  :: file_id     ! HDF file id's
      integer,        intent(out) :: nstate      ! number of statepoints on file

      logical :: ifxst
      integer :: ierror
      character(len=10) :: group_name

      nstate=0
      group_name='STATE_0001'

      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (.not. ifxst) then
        write (*,'(3a)') 'ERROR: Group STATE_0001 does not exist on file, is this a VERA output file?'
        return
      endif

      nstate=1    ! first state found successfully

      do
        write (group_name(7:10),'(i4.4)') nstate+1
        call h5lexists_f(file_id, group_name, ifxst, ierror)
        if (.not.ifxst) exit
        nstate=nstate+1
      enddo

      return
      end subroutine countstates

!=======================================================================
