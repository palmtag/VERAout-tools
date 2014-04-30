!=======================================================================
!
!  Program to read MPACT HDF output file and print summary
!
!  Copyright (c) 2014 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2014/04/28 - draft 1
!
!-----------------------------------------------------------------------
      program mpactread
      use  hdf5
      use  mod_hdftools
      implicit none

      character(len=80)  :: inputfile
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i, j, k, n
      integer            :: ierror
      integer            :: itype
      integer            :: ndim         ! temp variable
      integer            :: idim(10)     ! temp variable
      integer            :: nstate       ! statepoint number

      logical            :: ifxst
      logical            :: ifdebug=.false. ! debug flag

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=10)  :: group_name

      integer(hid_t)     :: file_id

! input data

      integer  :: isym                      ! core symmetry
      integer  :: iver                      ! output file version number
      integer  :: naxial                    ! number of axial levels
      integer  :: icore, jcore              ! core size
      integer, allocatable :: mapcore(:,:)  ! core map

      integer  :: kd                        ! number of axial levels in pin maps
      integer  :: nassm                     ! number of assemblies in pin maps
      integer  :: npin                      ! pin data

      real(8)  :: xkeff                     ! eigenvalue
      real(8)  :: pinmax                    ! peak 3D pin

      real(8), allocatable :: axial(:)      ! axial elevations
      real(8), allocatable :: pin2 (:,:,:)  ! 2d collapsed pin powers
      real(8), allocatable :: powertemp(:,:,:,:)
      real(8), allocatable :: power(:,:,:,:)

      character(len=2), allocatable :: xlabel(:)  ! assembly map labels
      character(len=2), allocatable :: ylabel(:)  ! assembly map labels

      integer, parameter :: maxstate=200
      real(8)  :: state_xkeff(maxstate)
      real(8)  :: state_pinmax(maxstate)

! command line flags

      logical  :: if3d   =.false.           ! turn on 3D edits
      logical  :: if2d   =.false.           ! turn on 2D edits
      logical  :: if2da  =.false.           ! turn on 2D assembly edits
      logical  :: if1d   =.false.           ! turn on 1D edits

!  initialize

      inputfile=' '
      state_xkeff(:)=0.0d0
      state_pinmax(:)=0.0d0

!----------------------------------------------------------------------
!  Read in arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  mpactread.exe [hdf5_file] {1D/2D/2DA/3D}'
        stop
      endif

! parse command line arguments

      do i=1, iargs
        call get_command_argument(i,carg)
        if     (carg.eq.'3D' .or. carg.eq.'3d') then
          if3d=.true.
        elseif (carg.eq.'2D' .or. carg.eq.'2d') then
          if2d=.true.
        elseif (carg.eq.'2DA' .or. carg.eq.'2da') then
          if2da=.true.
        elseif (carg.eq.'1D' .or. carg.eq.'1d') then
          if1d=.true.
        else
          inputfile=carg
        endif
      enddo

      if (.not.if1d)  write (*,*) 'no 1D edits requested on command line'
      if (.not.if2d)  write (*,*) 'no 2D edits requested on command line'
      if (.not.if2da) write (*,*) 'no 2DA edits requested on command line'
      if (.not.if3d)  write (*,*) 'no 3D edits requested on command line'

!--- initialize fortran interface

      call h5open_f(ierror)

!--- open HDF file

      write (*,'(2a)') 'reading h5 file: ', trim(inputfile)
      inquire(file=inputfile, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(inputfile),' does not exist'
        stop
      endif
      call h5fopen_f (inputfile, H5F_ACC_RDWR_F, file_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop
      endif

!--------------------------------------------------------------------------------
! Read top level HDF data
!--------------------------------------------------------------------------------
!  DATASET "axial_mesh" {
!  DATASET "core_map" {
!  DATASET "core_sym" {
!  DATASET "pin_volumes" {
!  DATASET "version" {

!? DATASET "keff" {
!? DATASET "pin_powers" {

!--------------------------------------------------------------------------------
! Read Eigenvalue and other scalers
!--------------------------------------------------------------------------------

      idim(:)=0
      naxial=0          ! number of axial edit bounds
      icore=0           ! number of assemblies across
      jcore=0           ! number of assemblies across

      xkeff=-100.0d0
      isym=-100
      iver=-100

!*** top level eigenvalue should not be here, it is a mistake
      dataset='keff'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call read_double(file_id, dataset, xkeff)  ! *** should not be here
        write (*,*) 'warning: found top level keff = ', xkeff
      endif

      dataset='core_sym'
      call read_integer(file_id, dataset, isym)

      dataset='version'
      call read_integer(file_id, dataset, iver)

      if (ifdebug) write (*,*) 'debug: isym  = ', isym
      if (ifdebug) write (*,*) 'debug: iver  = ', iver

      dataset='core_map'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.2) stop 'invalid dimensions in core_map'
      icore=idim(1)
      jcore=idim(2)
      if (ifdebug) write (*,*) 'debug: icore = ', icore
      if (ifdebug) write (*,*) 'debug: jcore = ', jcore
      allocate (mapcore(icore,jcore))
      call read_integer2d(file_id, dataset, icore, jcore, mapcore)

      write (*,'(/,a)') ' Core Map:'
      do j=1, jcore
         write (*,'(2x,20i3)') (mapcore(i,j),i=1,icore)
      enddo

!--- axial mesh

      dataset='axial_mesh'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.1) stop 'invalid dimensions in axial_mesh'
      naxial=idim(1)
      if (ifdebug) write (*,*) 'debug: naxial = ', naxial
      allocate (axial(naxial))
      call read_double1d(file_id, dataset, naxial, n, axial)   ! n will be set to naxial
      if (n.ne.naxial) stop 'error reading axial'

      write (*,'(/,a)') ' Axial Edit Boundaries'
      do j=1, naxial
        write (*,'(i5,20f12.4)') j, axial(j)
      enddo
      if (naxial.eq.0) then
        write (*,*) '*** warning: no axial edits found ****'
      endif

!--- temp define edit labels

  !  xlabel  R P N M L K J H G  F  E  D  C  B  A
  !  ylabel  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

      allocate (xlabel(icore))
      allocate (ylabel(jcore))

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

!***  dataset='pin_volumes'      ****

!--------------------------------------------------------------------------------
! Read statepoint data
!--------------------------------------------------------------------------------

! 40 format (/,'--------------------------',&
!            /,'  Reading statepoint ', i0, &
!            /,'--------------------------')

      dataset='/STATE_n/'
      nstate=0
      do
        nstate=nstate+1

        write (group_name,'(a,i0,a)') '/STATE_', nstate, '/'   ! *** needs to be i4.4 when mpact is fixed
        if (ifdebug) write (*,*) 'debug: state= ', group_name

!--- check if statepoint exists

        call h5lexists_f(file_id, group_name, ifxst, ierror)
        if (.not.ifxst) then
          if (ifdebug) write (*,*) 'dataset not found - exiting'
          goto 800
        endif

!!      write (*,40) nstate   ! statepoint is already printed in reading dataset edits

!--- eigenvalue

        dataset=trim(group_name)//'keff'
        call read_double(file_id, dataset, xkeff)

        if (ifdebug) write (*,*) 'debug: keff = ', xkeff

!--------------------------------------------------------------------------------
! Read pin powers
!--------------------------------------------------------------------------------

!--- allocate arrays if first statepoint

        if (nstate.eq.1) then
          dataset=trim(group_name)//'pin_powers'
          call h5info(file_id, dataset, itype, ndim, idim)

          nassm=idim(1)    ! order (nassm, kd, npin, npin)
          kd   =idim(2)
          npin =idim(4)
          if (ifdebug) then
            write (*,*) 'debug pin: ndim  = ', ndim
            write (*,*) 'debug pin: nassm = ', nassm
            write (*,*) 'debug pin: kd    = ', kd
            write (*,*) 'debug pin: npin  = ', npin
          endif

          if (ndim.ne.4) stop 'invalid dimensions in pin data'
          if (nassm.eq.0) stop 'invalid number of assemblies in pin data'
          if (npin .eq.0) stop 'invalid number of pins in pin data'
          if (idim(3).ne.idim(4)) stop 'invalid npin in pin data'

          if (naxial.ne.kd) then
            stop 'mismatch found between axial levels in input and number of levels of output'
          endif

          if (ifdebug) write (*,*) 'debug: allocating pin arrays'

          allocate (powertemp(nassm, kd, npin, npin))
          allocate (power    (npin, npin, kd, nassm))
          allocate (pin2     (npin, npin, nassm))
        endif

!--- array order was modified 12/6/2013 so that hdfview and h5dump look correct

        dataset=trim(group_name)//'pin_powers'
        call read_double4d(file_id, dataset, nassm, kd, npin, npin, powertemp)

!--- move power to old order

        do n=1, nassm
          do k=1, kd
            do j=1, npin
              do i=1, npin
                power(i,j,k,n)=powertemp(n,k,j,i)
              enddo
            enddo
          enddo
        enddo

!--- collapse

        write (*,'(a, f12.7)') ' keff =', xkeff

        call collapse(npin, kd, nassm, power, axial, pin2, pinmax)

!--- save statepoint values

        state_xkeff(nstate)=xkeff
        state_pinmax(nstate)=pinmax

!--- print maps

        if (if3d) call print_pin_map(npin, kd, power, nassm)

        if (if2d) call print_pin_map(npin, 1,  pin2, nassm)

        if (if2da) call print2d_assm_map(npin, nassm, pin2, icore, jcore, mapcore, xlabel, ylabel)

!!      if (if1d) call print1d(npin, kd, nassm, power, pin2, axial, icore, jcore, mapcore)

        write (*,*)

      enddo   ! end of statepoint loop

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------
  800 continue

      call h5fclose_f(file_id, ierror)

! print summary

      nstate=nstate-1   ! decrease due to statepoint check
      write (*,110)
      do n=1, nstate
        write (*,120) n, state_xkeff(n), state_pinmax(n)
      enddo
  110 format (/,'==================================',&
              /,'       Statepoint Summary', &
              /,'==================================',&
              /,'   N   eigenvalue   pinmax')

  120 format (i4, f12.5, f10.4)

! deallocate memory and stop

      if (allocated(pin2)) then  ! protect from missing statepoints
        deallocate (pin2)
        deallocate (power)
        deallocate (powertemp)
      endif

      deallocate (mapcore)
      deallocate (xlabel,ylabel)
      deallocate (axial)

      write (*,'(/,a)') 'done'

      end program
!=======================================================================
!
!  Collapse 3D edits to 2D and 1D
!
!=======================================================================
      subroutine collapse(npin, kd, nassm, power, axial, pin2, pinmax)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)
      real(8), intent(out):: pin2(npin,npin,nassm)
      real(8), intent(out):: pinmax    ! peak 3D pin

!--- local

      integer :: i, j, k
      integer :: ia
      integer :: kmin(4)
      integer :: kmax(4)
      real(8) :: pp
      real(8) :: zave
      real(8) :: cmin, cmax

      pin2(:,:,:)=0.0d0

!--- collapse pin powers to 2D

      cmax=0.0d0         ! core max
      cmin=1.0d20        ! core min
      kmax(:)=0          ! core max i, j, k, n
      kmin(:)=0          ! core min i, j, k, n

      do ia=1, nassm     ! loop over assemblies
        do j=1, npin
          do i=1, npin
            zave=0.0d0   ! axial height with power
            do k=1, kd   ! loop over axial levels
              pp=power(i,j,k,ia)
              if (pp.gt.0.0d0) then
                pin2(i,j,ia)=pin2(i,j,ia)+pp*axial(k)
                zave=zave+axial(k)
                if (pp.lt.cmin) then
                  cmin=pp
                  kmin(1)=i
                  kmin(2)=j
                  kmin(3)=k
                  kmin(4)=ia
                endif
                if (pp.gt.cmax) then
                  cmax=pp
                  kmax(1)=i
                  kmax(2)=j
                  kmax(3)=k
                  kmax(4)=ia
                endif
              endif
            enddo
            if (zave.gt.0.0d0) pin2(i,j,ia)=pin2(i,j,ia)/zave
          enddo
        enddo
      enddo      ! ia

      write (*,180) 'max', cmax, kmax(:)
      write (*,180) 'min', cmin, kmin(:)

      pinmax=cmax   ! save peak 3D pin

  180 format (' 3D ',a,' pin in core =', f10.4,' at (i,j,k,na)', 4i4)

      return
      end subroutine collapse

!=======================================================================
!
!  Print pin powers by assembly using pretty output
!
!=======================================================================
      subroutine print_pin_map(npin, kd, power, nassm)
!!    use mod_input, only : nassm, assmmap
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)

      integer            :: i, j, m
      integer            :: ia, klev, kpin
      real(8)            :: pp, pave, pmin, pmax
      character(len=150) :: line

      if (npin.gt.20) then
        write (*,*) '**** too many pins across to print nice maps ****'
        return
      endif

!**** To-do: read labels from input and use in edits

      do ia=1, nassm        ! loop over assemblies
        do klev=kd, 1, -1   ! loop over axial levels
!!        write (*,'(/,1x,a,i3,2a)') 'Assembly ', ia,'    type ', trim(assmmap(nassm))   ! no types to mpact
          write (*,'(/,1x,a,i3,2a)') 'Assembly ', ia
          if (kd.gt.1) then
             write (*,'(  1x,a,i3)') 'Level    ', klev
          else
             write (*,*) '2D collapse'
          endif
          pave=0.0
          kpin=0
          pmin=1.0d20
          pmax=0.0d0
          do i=1, npin
            line=' '
            m=0
            do j=1, npin
              pp=power(i,j,klev,ia)
              if (pp.eq.0.0) then
                line(m+1:m+7)='   --- '
              else
                write (line(m+1:m+7),'(f7.4)') pp
                pave=pave+pp
                kpin=kpin+1
                pmin=min(pmin,pp)
                pmax=max(pmax,pp)
              endif
              m=m+7
            enddo
            write (*,'(1x,a)') line(1:m)
          enddo
          if (kpin.gt.0) pave=pave/real(kpin)
          write (*,210) kpin, pave, pmin, pmax

        enddo    ! klev
      enddo      ! ia
  210 format (5x,'number of hot pins',i4,'    average=',f7.4,'   min=', f7.4,'   max=', f7.4)

      return
      end subroutine print_pin_map
!=======================================================================
!
!  Subroutine to print 2D Assembly Maps
!
!=======================================================================
      subroutine print2d_assm_map(npin, nassm, pow2, icore, jcore, mapcore, xlabel, ylabel)
!!    use mod_input, only : icoresize, mapcore, xlabel, ylabel
      implicit none
      integer, intent(in) :: npin, nassm
      integer, intent(in) :: icore, jcore
      integer, intent(in) :: mapcore(icore,jcore)
      real(8), intent(in) :: pow2(npin,npin,nassm)
      character(len=*), intent(in) :: xlabel(icore)
      character(len=*), intent(in) :: ylabel(jcore)

      integer          :: i, j, nn
      integer          :: nnsave
      integer          :: ia, ja, na
      real(8)          :: pmin, pmax
      real(8)          :: pp
      real(8)          :: passm(nassm)  ! automatic

      write (*,*)

!--- calculate 2D assembly averages

      passm(:)=0.0d0

      nnsave=0
      do na=1, nassm
        nn=0
        do j=1, npin
          do i=1, npin
            pp=pow2(i,j,na)
            if (pp.gt.0.0d0) then
              nn=nn+1
              passm(na)=passm(na)+pp
            endif
          enddo
        enddo
        if (nn.gt.0) passm(na)=passm(na)/dble(nn)
        if (nnsave.eq.0) nnsave=nn    ! save first time
        if (nn.ne.nnsave) then
           write (*,*) 'ERROR: Assemblies have different pin counts', nn, nnsave
        endif
      enddo

!--- calculate average - use mapcore in case this is qtr-core

      nn=0
      pp=0.0d0
      pmin=passm(1)
      pmax=passm(1)

      do ja=1, jcore
        do ia=1, icore
          na=mapcore(ia,ja)
          if (na.gt.0) then
            pp=pp+passm(na)
            nn=nn+1
            pmax=max(pmax,passm(na))
            pmin=min(pmin,passm(na))
          endif
        enddo    ! ia
      enddo      ! ja

      if (nn.gt.0) pp=pp/dble(nn)
      write (*,*) 'number of assemblies in full-core', nn
      write (*,130) 'average', pp
      if (abs(pp-1.0d0).gt.0.0001) write (0,*) '***** check normalization *****'
      write (*,130) 'maximum', pmax
      write (*,130) 'minimum', pmin
  130 format (1x,a,' assembly power  ', f8.4)

!--- print map

      write (*,'("  **  ",50(4x,a2,2x))') (xlabel(i),i=1,icore)   ! labels are 2

      do j=1, jcore
        write (*,'(2x,a2,"- ")',advance='no') ylabel(j)
        do i=1, icore
           if (mapcore(i,j).eq.0) then
             write (*,'(8x)',advance='no')
           else
             write (*,'(f8.4)',advance='no') passm(mapcore(i,j))
           endif
        enddo
        write (*,*)
      enddo

      return
      end subroutine print2d_assm_map
