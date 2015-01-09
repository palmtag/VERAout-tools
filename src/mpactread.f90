   program mpactread
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
!  2014/05/01 - update for VERAOUT version 1
!  2014/08/01 - change HDF open statement to read-only
!             - read exposure values
!
!-----------------------------------------------------------------------
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
      integer            :: nstate=0     ! statepoint number

      logical            :: ifxst
      logical            :: ifdebug=.false. ! debug flag
      logical            :: ifoldstate=.false.

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=12)  :: group_name

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

      real(8)  :: xtemp                     ! temp value
      real(8)  :: xkeff                     ! eigenvalue
      real(8)  :: xexpo                     ! exposure
      real(8)  :: pinmax                    ! peak 3D pin
      real(8)  :: rated_power               ! rated flow
      real(8)  :: rated_flow                ! rated power

      real(8), allocatable :: axial(:)      ! axial elevations
      real(8), allocatable :: pin2 (:,:,:)  ! 2d collapsed pin powers
      real(8), allocatable :: powertemp(:,:,:,:)
      real(8), allocatable :: power(:,:,:,:)

      character(len=80)  :: title           ! Problem title
      character(len=2), allocatable :: xlabel(:)  ! assembly map labels
      character(len=2), allocatable :: ylabel(:)  ! assembly map labels

      integer, parameter :: maxstate=200
      real(8)  :: state_xkeff(maxstate)
      real(8)  :: state_xexpo(maxstate)  ! exposure
      real(8)  :: state_pinmax(maxstate)

! command line flags

      logical  :: if3d   =.false.           ! turn on 3D edits
      logical  :: if2d   =.false.           ! turn on 2D edits
      logical  :: if2da  =.false.           ! turn on 2D assembly edits
      logical  :: if1d   =.false.           ! turn on 1D edits

!  initialize

      inputfile=' '
      state_xkeff(:)=0.0d0
      state_xexpo(:)=0.0d0
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
      call h5fopen_f (inputfile, H5F_ACC_RDONLY_F, file_id, ierror)   ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop
      endif

!--------------------------------------------------------------------------------
! Read top level HDF data
!--------------------------------------------------------------------------------

      xkeff=-100.0d0
      xexpo=-100.0d0
      iver=-100

!*** top level eigenvalue should not be here, it is a mistake
      dataset='keff'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (ifxst) then
        call hdf5_read_double(file_id, dataset, xkeff)  ! *** should not be here
        write (*,*) 'warning: found top level keff = ', xkeff
      endif

! version

      dataset='veraout_version'
      call hdf5_read_integer(file_id, dataset, iver)

! title

      title=' '

      dataset='title'
      call read_string(file_id, dataset, title)

      write (*,*)
      write (*,'(2a)')     ' Title: ', trim(title)
      write (*,'(a,i2)') ' VERAOUT version ', iver
      write (*,*)

!-------------------
!  Read CORE group
!-------------------

      isym=-100
      idim(:)=0
      naxial=0          ! number of axial edit bounds
      icore=0           ! number of assemblies across
      jcore=0           ! number of assemblies across

! check of CORE group exists (group may not exist on old MPACT files and all data will be in root)

      group_name='/CORE/'
      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (.not.ifxst) then
        group_name=' '
        write (*,*) 'CORE group not found - is this an old MPACT file?'
      endif

!! no need to open group on read???

      dataset=trim(group_name)//'core_sym'
      call hdf5_read_integer(file_id, dataset, isym)

!--- rated power and flow

      rated_power=0.0d0   ! default
      rated_flow =0.0d0   ! default

      dataset=trim(group_name)//'rated_power'
      call hdf5_read_double(file_id, dataset, rated_power)

      dataset=trim(group_name)//'rated_flow'
      call hdf5_read_double(file_id, dataset, rated_flow)

      write (*,*)
      write (*,'(a,i2)')    ' core symmetry  isym  = ', isym
      write (*,'(a,f12.4)') ' rated power ', rated_power
      write (*,'(a,f12.4)') ' rated flow  ', rated_flow

!--- core map

      dataset=trim(group_name)//'core_map'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.2) stop 'invalid dimensions in core_map'
      icore=idim(1)
      jcore=idim(2)
      if (ifdebug) write (*,*) 'debug: icore = ', icore
      if (ifdebug) write (*,*) 'debug: jcore = ', jcore
      allocate (mapcore(icore,jcore))
      call hdf5_read_integer(file_id, dataset, icore, jcore, mapcore)

      write (*,'(/,a)') ' Core Map:'
      do j=1, jcore
         write (*,'(2x,20i3)') (mapcore(i,j),i=1,icore)
      enddo

!--- axial mesh

      dataset=trim(group_name)//'axial_mesh'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.1) stop 'invalid dimensions in axial_mesh'
      naxial=idim(1)   ! actually naxial+1
      if (ifdebug) write (*,*) 'debug: naxial+1 = ', naxial
      allocate (axial(naxial))
      call hdf5_read_double(file_id, dataset, naxial, n, axial)   ! n will be set to naxial
      if (n.ne.naxial) then
         write (*,*) 'found    ', n
         write (*,*) 'expecting', naxial+1
         stop 'error reading axial'
      endif
      naxial=naxial-1    ! decrease to match number of levels instead of boundaries

      if (naxial.eq.0) then   ! fixup OLD MPACT
        write (*,*) 'ERROR in axial edits - is this an old broken MPACT file?'
        naxial=1
        axial(1)=0.0d0
        axial(2)=1.0d0
      endif

      write (*,'(/,a)') ' Axial Edit Boundaries'
      do j=1, naxial+1
        write (*,'(i5,20f12.4)') j-1, axial(j)
      enddo
      if (naxial.eq.0) then
        write (*,*) '*** warning: no axial edits found ****'
      endif

      write (*,'(a,f12.4)') ' total core height = ', axial(naxial+1)-axial(1)
      write (*,*)

      do j=1, naxial
        axial(j)=axial(j+1)-axial(j)   ! convert to deltas
      enddo
      axial(naxial+1)=0.0d0

      if (ifdebug) then
        write (*,*) 'axial deltas'
        do j=1, naxial+1
          write (*,'(i5,20f12.4)') j, axial(j)
        enddo
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

!--- close group (if opened)

!!    if (group_name.ne.' ') then
!!      call h5gclose_f (group_id, ierror)
!!      if (ierror.ne.0) stop 'error closing group'
!!    endif


!---------------------
!  Read STATE groups
!---------------------

!--- check if this file uses an old mpact STATE format

      group_name='/STATE_1/'
      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (ifxst) ifoldstate=.true.

! 40 format (/,'--------------------------',&
!            /,'  Reading statepoint ', i0, &
!            /,'--------------------------')

      nstate=0
      do
        nstate=nstate+1

        if (ifoldstate) then
          write (group_name,'(a,i0,a)') '/STATE_', nstate, '/'
        else
          write (group_name,'(a,i4.4,a)') '/STATE_', nstate, '/'
        endif
        if (ifdebug) write (*,*) 'debug: state= ', group_name

!--- check if statepoint exists

        call h5lexists_f(file_id, group_name, ifxst, ierror)
        if (.not.ifxst) then
          if (ifdebug) write (*,*) 'dataset not found - exiting'
          goto 800
        endif

!!      write (*,40) nstate   ! statepoint is already printed in reading dataset edits

!--- eigenvalue

! read exposure value at this statepoint

!alt    dataset=trim(group_name)//'exposure_GWDMT'
        dataset=trim(group_name)//'exposure'
        call hdf5_read_double(file_id, dataset, xexpo)

        dataset=trim(group_name)//'keff'
        call hdf5_read_double(file_id, dataset, xkeff)

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
        call hdf5_read_double(file_id, dataset, nassm, kd, npin, npin, powertemp)

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

        write (*,'(a, f10.4,a)') ' exposure =', xexpo,' GWD/MT'
        write (*,'(a, f12.7)') ' keff =', xkeff

        call collapse(npin, kd, nassm, power, axial, pin2, pinmax, icore, jcore, mapcore)

!--- save statepoint values

        if (nstate.gt.maxstate) stop 'maxstate exceeded - increase and recompile'

        state_xkeff(nstate)=xkeff
        state_xexpo(nstate)=xexpo
        state_pinmax(nstate)=pinmax

!--- print maps

        call stat3d('pin_powers', npin,  kd, nassm, axial, power, xtemp)

        if (if3d) call print_3D_pin_map('3D pin_powers', npin, kd, nassm, power)

        if (if2d) call print_3D_pin_map('2D pin_powers', npin, 1,  nassm, pin2)

        if (if2da) call print2d_assm_map(npin, nassm, pin2, icore, jcore, mapcore, xlabel, ylabel)

        if (if1d) then
          call print1d('pin_powers', npin, kd, nassm, power, axial)
        endif

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
        write (*,120) n, state_xexpo(n), state_xkeff(n), state_pinmax(n)
      enddo
  110 format (/,'==================================',&
              /,'       Statepoint Summary', &
              /,'==================================',&
              /,'   N   exposure  eigenvalue   pinmax')

  120 format (i4, f10.4, f12.5, f10.4)

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
      subroutine collapse(npin, kd, nassm, power, axial, pin2, pinmax, icore, jcore, mapcore)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      integer, intent(in) :: icore, jcore
      integer, intent(in) :: mapcore(icore,jcore)
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)
      real(8), intent(out):: pin2(npin,npin,nassm)
      real(8), intent(out):: pinmax    ! peak 3D pin

!--- local

      integer :: i, j, k
      integer :: ia, ja, na
      integer :: k3min(4)
      integer :: k3max(4)
      integer :: k2min(3)
      integer :: k2max(3)
      real(8) :: pp
      real(8) :: zave, zrod
      real(8) :: c3min, c3max
      real(8) :: c2min, c2max
      real(8) :: cave, clen

      pin2(:,:,:)=0.0d0

      c3max=0.0d0        ! 3D core max
      c3min=1.0d20       ! 3D core min
      c2max=0.0d0        ! 2D core max
      c2min=1.0d20       ! 2D core min
      cave=0.0d0         ! core average
      clen=0.0d0         ! core fuel length
      k3max(:)=0         ! 3D core max i, j, k, n
      k3min(:)=0         ! 3D core min i, j, k, n
      k2max(:)=0         ! 2D core max i, j, n
      k2min(:)=0         ! 2D core min i, j, n

      write (*,'(/,1x,a)') 'Collapsing 3D edits to 2D'

!--- check 3D normalization - loop over mapcore in case symmetry was used


      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        do j=1, npin
          do i=1, npin
            do k=1, kd   ! loop over axial levels
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                cave=cave+axial(k)*pp
                clen=clen+axial(k)
                if (pp.lt.c3min) then
                  c3min=pp
                  k3min(1)=i
                  k3min(2)=j
                  k3min(3)=k
                  k3min(4)=na
                endif
                if (pp.gt.c3max) then
                  c3max=pp
                  k3max(1)=i
                  k3max(2)=j
                  k3max(3)=k
                  k3max(4)=na
                endif
              endif
            enddo
          enddo
        enddo

      enddo      ! ia
      enddo      ! ja

      if (clen.gt.0.0d0) cave=cave/clen

      if (abs(cave-1.0d0).gt.1.0d-6) then
         write (*,*) 'WARNING: pin powers do not appear to be normalized correctly'
         write (*,*) 'total core length ', clen
         write (*,*) 'core average power of linear power ', cave
      endif

      write (*,180) 'max', c3max, k3max(:)
      write (*,180) 'min', c3min, k3min(:)

      pinmax=c3max   ! save peak 3D pin

  180 format (' 3D ',a,' pin in core =', f10.4,' at (i,j,k,na)', 4i4)
  190 format (' 2D ',a,' pin in core =', f10.4,' at (i,j,na)  ', 2i4,4x,i4)

!--- collapse pin powers to 2D

      do na=1, nassm     ! loop over assemblies
        do j=1, npin
          do i=1, npin
            zrod=0.0d0
            zave=0.0d0   ! 2D axial height with power
            do k=1, kd   ! loop over axial levels
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                zrod=zrod+axial(k)*pp
                zave=zave+axial(k)
              endif
            enddo
            if (zave.gt.0.0d0) then
              pin2(i,j,na)=zrod/zave
              pp=pin2(i,j,na)
              if (pp.lt.c2min) then
                c2min=pp
                k2min(1)=i
                k2min(2)=j
                k2min(3)=na
              endif
              if (pp.gt.c2max) then
                c2max=pp
                k2max(1)=i
                k2max(2)=j
                k2max(3)=na
              endif
            endif
          enddo
        enddo
      enddo      ! na

      write (*,190) 'max', c2max, k2max(:)
      write (*,190) 'min', c2min, k2min(:)

      return
      end subroutine collapse

!=======================================================================
!
!  Subroutine to print 2D Assembly Maps
!
!=======================================================================
      subroutine print2d_assm_map(npin, nassm, pow2, icore, jcore, mapcore, xlabel, ylabel)
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

      logical  :: ifbw  ! flag if assemblies have different number of pins

      ifbw=.false.

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
        if (nn.ne.nnsave) ifbw=.true.
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
      if (ifbw) then
        write (*,*) 'WARNING: average of assemblies cannot be calculated because assemblies have different number of rods'
      else
        write (*,130) 'average', pp
        if (abs(pp-1.0d0).gt.0.0001) write (0,*) '***** check normalization *****'
      endif
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
