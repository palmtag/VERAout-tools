   program veradiff2
   use hdf5
   implicit none
!=======================================================================
!
!  Program to compare two VERA HDF output file and print summary
!
!  *** This version allows the comparison of two files with different numbers of exposure points ***
!
!  *** Use veradiff_testing for automated testing applications ****
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
!  2017/05/01 - add ability to test files with different exposure points
!             - created new version "veradiff_testing" for automated testing applications
!
!--------------------------------------------------------------------------------

      character(len=250) :: fname1, fname2   ! HDF file names
      integer            :: iargs            ! number of command line arguments

      integer(hid_t)     :: file_id1, file_id2  ! HDF file id's
      integer            :: ierror
      integer            :: i, n, n2
      integer            :: nstate           ! statepoint number
      integer            :: nstate2          ! statepoint number

      integer, parameter :: maxstate=200

      real(8) :: xkeff1(maxstate)
      real(8) :: xkeff2(maxstate)
      real(8) :: xexpo1(maxstate)
      real(8) :: xexpo2(maxstate)
      real(8) :: xefpd1(maxstate)
      real(8) :: xefpd2(maxstate)
      real(8) :: xboron1(maxstate)
      real(8) :: xboron2(maxstate)
      real(8) :: dk, db
      real(8) :: xmax(maxstate)
      real(8) :: xrms(maxstate)

      logical            :: ifxst

      ierror=0
      fname1=' '
      fname2=' '

      nstate=0   ! number of statepoints on file 1
      nstate2=0  ! number of statepoints on file 2

      xexpo1=0.0d0
      xexpo2=0.0d0

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.2) then
        write (*,*) 'usage:  veradiff.exe [file1] [file2]'
        goto 900
      endif

      call get_command_argument(1,fname1)
      call get_command_argument(2,fname2)

      write (*,'(/,1x,2a)') 'reading h5 file: ', trim(fname1)
      inquire(file=fname1, exist=ifxst)
      if (.not.ifxst) then
        write (*,*) 'ERROR: input file ',trim(fname1),' does not exist'
        goto 900
      endif

      write (*,'(1x,2a)') 'reading h5 file: ', trim(fname2)
      inquire(file=fname2, exist=ifxst)
      if (.not.ifxst) then
        write (*,*) 'ERROR: input file ',trim(fname2),' does not exist'
        goto 900
      endif

!--------------------------------------------------------------------------------
! Read HDF files
!--------------------------------------------------------------------------------

!--- initialize fortran interface

      call h5open_f(ierror)      ! NOTE: THIS IS REQUIRED

!--- open file1

      call h5fopen_f (fname1, H5F_ACC_RDWR_F, file_id1, ierror)
      if (ierror<0) then
        write (*,'(3a)') 'ERROR: H5 input file ',trim(fname1),' could not be opened'
        goto 900
      endif

!--- open file2

      call h5fopen_f (fname2, H5F_ACC_RDWR_F, file_id2, ierror)
      if (ierror<0) then
        write (*,'(3a)') 'ERROR: H5 input file ',trim(fname2),' could not be opened'
        goto 900
      endif

!--- read number of statepoints and eigenvalues

      call countstates(file_id1, maxstate, nstate, xexpo1, xefpd1, xkeff1, xboron1)

      call countstates(file_id2, maxstate, nstate2, xexpo2, xefpd2, xkeff2, xboron2)

!---- calculate pin difference values (experimental)

      do n=1, nstate
        n2=0
        do i=1, nstate2
          if (abs(xexpo1(n)-xexpo2(i)).lt.0.0001d0) n2=i
        enddo
        if (n2.eq.0) cycle  ! no match found

        call pin_compare(file_id1, file_id2, n, n2, xmax(n), xrms(n))

      enddo

!--------------------------------------------------------------------------------
! Print Results
!--------------------------------------------------------------------------------

      write (*,*)
      write (*,'(2a)') ' file1: ', trim(fname1)
      write (*,'(2a)') ' file2: ', trim(fname2)

      write (*,122)

      do n=1, nstate
        n2=0
        do i=1, nstate2 
          if (abs(xexpo1(n)-xexpo2(i)).lt.0.0001d0) n2=i
        enddo
        if (n2.eq.0) cycle  ! no match found

        dk=(xkeff2(n2)-xkeff1(n))*1.0d5
        db= xboron2(n2)-xboron1(n)

        write (*,55) n, xexpo1(n), xefpd1(n), xkeff1(n), xkeff2(n2), dk, xboron1(n), xboron2(n2), db, xmax(n), xrms(n)

      enddo

  122 format ( &
            '  ==============================================================================', &
          /,'                         Statepoint Summary (file2-file1)', &
          /,'  ==============================================================================', &
          /,'     N     exp     EFPD    keff1     keff2      pcm     boron1   boron2     dbor   pinmax   pinrms')
   55 format (2x,i4, f9.4, f8.2, 2f10.6, f8.1, 2x, 3f9.2, 2f9.2)

!--- close HDF files

      call h5fclose_f(file_id1, ierror)
      call h5fclose_f(file_id2, ierror)

  900 continue

      end program

!=======================================================================
!
!  Subroutine to read two pin_power datasets and calculate RMS
!
!=======================================================================
      subroutine pin_compare(file_id1, file_id2, nstate1, nstate2, xmax, xrms)
      use hdf5
      use mod_hdftools, only : h5info, hdf5_read_double
      implicit none
      integer(hid_t),   intent(in)  :: file_id1, file_id2
      integer,          intent(in)  :: nstate1, nstate2
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

      character(len=24)  :: dataset1
      character(len=24)  :: dataset2

!--- dataset names

      dataset1='STATE_0001/pin_powers'
      dataset2='STATE_0001/pin_powers'

      write (dataset1(7:10),'(i4.4)') nstate1
      write (dataset2(7:10),'(i4.4)') nstate2

!d    write (*,*) 'debug: ', dataset1, nstate1
!d    write (*,*) 'debug: ', dataset2, nstate2

!--- find dimensions so we can allocate

      xmax=999999.0d0   ! set high in case reads fail
      ifbad=.false.

      call h5info(file_id1, dataset1, itype, ndim, idim)

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

      call hdf5_read_double(file_id1, dataset1, idim(1), idim(2), idim(3), idim(4), power1)
      call hdf5_read_double(file_id2, dataset2, idim(1), idim(2), idim(3), idim(4), power2)

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
 125  format (3x,'ave power',i1,5x, f10.4,'   (does not include volume weighting or qtr sym)')
 126  format (3x,'max power',i1,5x, f10.4,'   at (ip,jp,k,na) ', 4i4)
 130  format (3x,'rms difference ', f10.4,' %')
 131  format (3x,'max difference ', f10.4,' % at (ip,jp,k,na) ', 4i4)
 138  format (3x,'difference in max pin ', f10.4,' %')


      return
      end subroutine pin_compare
!=======================================================================
!
!  Subroutine to count number of statepoints on a file and read some scalars
!
!=======================================================================
      subroutine countstates(file_id, maxstate, nstate, xexpo, xefpd, xkeff, xboron)
      use mod_hdftools, only : hdf5_read_double
      use hdf5
      implicit none
      integer(hid_t), intent(in)  :: file_id     ! HDF file id's
      integer,        intent(in)  :: maxstate    ! max number of statepoints on file
      integer,        intent(out) :: nstate      ! number of statepoints on file
      real(8)                     :: xexpo(maxstate)
      real(8)                     :: xefpd(maxstate)
      real(8)                     :: xkeff(maxstate)
      real(8)                     :: xboron(maxstate)

      logical :: ifxst
      integer :: ierror
      character(len=80) :: dataset
      character(len=10) :: statename

      nstate=0
      statename='STATE_0001'

      do
        write (statename(7:10),'(i4.4)') nstate+1
        call h5lexists_f(file_id, statename, ifxst, ierror)
        if (.not.ifxst) exit
        nstate=nstate+1
        if (nstate.gt.maxstate) stop 'maxstate exceeded'

        dataset=trim(statename)//'/exposure'
        call hdf5_read_double(file_id, dataset, xexpo(nstate))

        dataset=trim(statename)//'/exposure_efpd'
        call hdf5_read_double(file_id, dataset, xefpd(nstate))

        dataset=trim(statename)//'/keff'
        call hdf5_read_double(file_id, dataset, xkeff(nstate))

        dataset=trim(statename)//'/boron'
        call hdf5_read_double(file_id, dataset, xboron(nstate))

      enddo

      return
      end subroutine countstates

!=======================================================================
