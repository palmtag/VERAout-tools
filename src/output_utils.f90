   subroutine stat3d (title, npin, kd, nassm, power, pave, pmax, pmax2d)
   use mod_coregeom, only : icore, jcore, mapcore, axial
   implicit none
!=======================================================================
!
!  Subroutines for writing 1D and 2D maps to output
!
!  Copyright (c) 2014-2017 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  Subroutine to print 3D statistics
!
!  General version that considers qtr-symmetry
!
!  2017/05/18 - added 2PIN edit
!
!=======================================================================
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(out):: pave   ! average
      real(8), intent(out):: pmax   ! max 3pin
      real(8), intent(out):: pmax2d ! max 2pin
      character(len=*), intent(in) :: title

!--- local

      integer  :: i, j, k
      integer  :: ia, ja, na, kpin, nexit
      real(8)  :: p2d, z2d
      real(8)  :: pp
      real(8)  :: zlen, zave      ! axial values
      real(8)  :: c3min, c3max    ! 3D values
      real(8)  :: exave           ! exit average
      integer  :: k2max(3)        ! 2D max locations
      integer  :: k3min(4)        ! 3D min locations
      integer  :: k3max(4)        ! 3D max locations

!--- calculate 3D statistics

      pave=0.0d0   ! output average
      pmax=0.0d0   ! output maximum 3pin
      pmax2d=0.0d0 ! output maximum 2pin

      kpin=0
      zave=0.0d0
      zlen=0.0d0
      c3min=1.0d20   ! core min/max
      c3max=0.0d0
      k3min(:)=0
      k3max(:)=0
      k2max(:)=0

      exave=0.0d0
      nexit=0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        do j=1, npin
          do i=1, npin
            p2d=0.0d0   ! reset for each rod
            z2d=0.0d0
            do k=1, kd
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                p2d=p2d+axial(k)*pp
                z2d=z2d+axial(k)
                zave=zave+axial(k)*pp
                zlen=zlen+axial(k)
                kpin=kpin+1
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
            if (z2d.gt.0.0d0) then
              p2d=p2d/z2d
              if (p2d.gt.pmax2d) then
                pmax2d=p2d
                k2max(1)=i
                k2max(2)=j
                k2max(3)=na
              endif
            endif
            pp=power(i,j,kd,na)   ! exit value
            exave=exave+pp
            nexit=nexit+1
          enddo
        enddo
      enddo    ! ia
      enddo    ! ja
      if (zlen.gt.0.0d0) then
        zave=zave/zlen
      else
        c3min=0.0d0    ! avoid overflow
      endif

      write (*,'(/,1x,2a)') '3D Statistics - ', trim(title)

      if (c3max.gt.0.0d0 .and. c3max.lt.0.001d0) then ! print exponent
        write (*,190) 'Min  ', c3min
        write (*,190) 'Max  ', c3max
        write (*,192) 'Max2d', pmax2d
        write (*,192) 'Ave  ', zave
      else
        write (*,180) 'Min  ', c3min, k3min(:)
        write (*,180) 'Max  ', c3max, k3max(:)
        write (*,181) 'Max2D', pmax2d, k2max(:)
        write (*,182) 'Ave  ', zave
      endif
      write (*,'(a,i10)  ')  '  Num  ', kpin

      write (*,*) 'exit average ', exave/dble(nexit), nexit

  180 format (2x, a,f10.4,'  at (i,j,k,na)', 4i4)
  181 format (2x, a,f10.4,'  at (i,j,  na)', 2i4,4x,i4)
  182 format (2x, a,f10.4)
  190 format (2x, a,1p,e14.5,'  at (i,j,k,na)', 0p, 4i4)
  192 format (2x, a,1p,e14.5)

      pave=zave    ! return max value
      pmax=c3max   ! return max value

!--- check for very small values

!     do k=1, kd
!       do na=1, nassm
!         do j=1, npin
!           do i=1, npin
!             pp=power(i,j,k,na)
!             if (pp.gt.0.0d0) then
!               if (pp.lt.zave*1.0d-4) then
!                  write (*,*) 'small value detected ', i, j, k, na, pp
!               endif
!             endif
!           enddo
!         enddo
!       enddo
!     enddo

      return
      end subroutine stat3d
!=======================================================================
!
!  Subroutine to print 3D statistics
!
!  Simple version - does straight average without considering qtr-symmetry
!
!=======================================================================
      subroutine stat3d_simple (title, npin, kd, nassm, axial, power, pmax)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(out):: pmax
      character(len=*), intent(in) :: title

!--- local

      integer  :: i, j, k
      integer  :: ia, kpin
      real(8)  :: pp
      real(8)  :: zlen, zave      ! axial values
      real(8)  :: c3min, c3max    ! 3D values
      integer  :: k3min(4)        ! 3D min locations
      integer  :: k3max(4)        ! 3D max locations

!--- calculate 3D statistics

      pmax=0.0d0   ! output maximum

      kpin=0
      zave=0.0d0
      zlen=0.0d0
      c3min=1.0d20
      c3max=0.0d0
      k3min(:)=0
      k3max(:)=0

      do k=1, kd
        do ia=1, nassm
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,ia)
              if (pp.gt.0.0d0) then
                zave=zave+axial(k)*pp
                zlen=zlen+axial(k)
                kpin=kpin+1
                if (pp.lt.c3min) then
                  c3min=pp
                  k3min(1)=i
                  k3min(2)=j
                  k3min(3)=k
                  k3min(4)=ia
                endif
                if (pp.gt.c3max) then
                  c3max=pp
                  k3max(1)=i
                  k3max(2)=j
                  k3max(3)=k
                  k3max(4)=ia
                endif

              endif
            enddo
          enddo
        enddo
      enddo
      if (zlen.gt.0.0d0) then
        zave=zave/zlen
      else
        c3min=0.0d0    ! avoid overflow
      endif

      write (*,'(/,1x,2a)') '3D Statistics - ', trim(title)

      if (c3max.gt.0.0d0 .and. c3max.lt.0.001d0) then ! print exponent
        write (*,190) 'Min ', c3min
        write (*,190) 'Max ', c3max
        write (*,192) 'Ave ', zave
      else
        write (*,180) 'Min ', c3min, k3min(:)
        write (*,180) 'Max ', c3max, k3max(:)
        write (*,182) 'Ave ', zave
      endif
      write (*,'(a,i10)  ')  '  Num ', kpin

  180 format (2x, a,f10.4,'  at (i,j,k,na)', 4i4)
  182 format (2x, a,f10.4)
  190 format (2x, a,1p,e14.5,'  at (i,j,k,na)', 0p, 4i4)
  192 format (2x, a,1p,e14.5)

      pmax=c3max   ! return max value

!--- check for very small values

!     do k=1, kd
!       do ia=1, nassm
!         do j=1, npin
!           do i=1, npin
!             pp=power(i,j,k,ia)
!             if (pp.gt.0.0d0) then
!               if (pp.lt.zave*1.0d-4) then
!                  write (*,*) 'small value detected ', i, j, k, ia, pp
!               endif
!             endif
!           enddo
!         enddo
!       enddo
!     enddo


      return
      end subroutine stat3d_simple

!=======================================================================
!
!  Subroutine to print 3D exit values by assembly using pretty output
!  (useful for looking at exit temperatures/densities)
!
!=======================================================================
      subroutine print_exit_map(title, npin, kd, nassm, power)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer            :: ia, klev

!--- print exit values (level kd)

      klev=kd

      do ia=1, nassm        ! loop over assemblies
        write (*,*)
        write (*,*) 'Exit values level ', kd
        call print_single_pin_map(title, npin, kd, nassm, ia, klev, power)
      enddo      ! ia

      return
      end subroutine print_exit_map
!=======================================================================
!
!  Subroutine to print 3D pin powers by assembly using pretty output
!
!=======================================================================
      subroutine print_3D_pin_map(title, npin, kd, nassm, power)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer            :: ia, klev

!--- start

      do ia=1, nassm        ! loop over assemblies
        do klev=kd, 1, -1   ! loop over axial levels
          call print_single_pin_map(title, npin, kd, nassm, ia, klev, power)
        enddo    ! klev
      enddo      ! ia

      return
      end subroutine print_3D_pin_map
!=======================================================================
!
!  Subroutine to print a single pin power map of an assembly using pretty output
!
!  Pass in entire 3D array (nassm,kd), but only print one assembly (ia,klev)
!
!=======================================================================
      subroutine print_single_pin_map(title, npin, kd, nassm, ia, klev, power)
      implicit none
      integer, intent(in) :: npin, kd, nassm, klev, ia
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer            :: i, j, m
      integer            :: kpin
      real(8)            :: pp, pave, pmin, pmax
      real(8)            :: zmax
      character(len=200) :: line
      character(len=8)   :: fmt

      if (npin.gt.30) then
        write (*,*) '**** too many pins across to print nice maps ****'
        return
      endif

!--- calculate maximum and set the format correctly

      pp=0.0d0
      zmax=0.0d0
      do j=1, npin
        do i=1, npin
          pp=power(i,j,klev,ia)
          zmax=max(zmax,pp)
        enddo
      enddo

      fmt='(f7.4)'
      if (zmax.gt.  10.0d0) fmt='(f7.3)'
      if (zmax.gt. 100.0d0) fmt='(f7.2)'
      if (zmax.gt.1000.0d0) fmt='(f7.1)'

!d    write (*,*) 'debug: zmax = ', zmax
!d    write (*,*) 'debug: fmt  = ', fmt

!--- write map

      write (*,'(/,1x,a)') trim(title)
!!    write (*,'(1x,a,i3,2a)') 'Assembly ', ia,'    type ', trim(assmmap(nassm))   ! no types to mpact
      if (kd.eq.1) then
        write (*,30)  ia
      else
        write (*,35)  ia, klev
      endif
      pave=0.0
      kpin=0
      pmin=2.0d20
      pmax=0.0d0
      do j=1, npin
        line=' '
        m=0
        do i=1, npin
          pp=power(i,j,klev,ia)
          if (pp.eq.0.0) then
            line(m+1:m+7)='   --- '
          else
            write (line(m+1:m+7),fmt) pp
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
      if (pmin.gt.1.0d20) pmin=0.0d0   ! protect if no data
      write (*,210) kpin, pave, pmin, pmax
!d    write (*,*) 'debug average = ', pave

! 210 format (5x,'number of hot pins',i4,'    average=',f7.4,'   min=', f7.4,'   max=', f7.4)
  210 format (5x,'number of hot pins',i4,'    average=',f12.4,'   min=', f12.4,'   max=', f12.4)

   30 format (1x,'Assembly ',i3,2a)
   35 format (1x,'Assembly ',i3,'  Level',i3)

      return
      end subroutine print_single_pin_map
!=======================================================================
!
!  Subroutine to print axial edits
!
!=======================================================================
      subroutine print1d(title, npin, kd, nassm, power)
      use mod_coregeom, only : icore, jcore, mapcore, axial
      implicit none
      character(len=*), intent(in) :: title
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)

      integer          :: i, j, k, nn
      integer          :: ia, ja
      integer          :: na
      real(8)          :: z1, z2
      real(8)          :: pp
      real(8)          :: pave
      real(8)          :: zave
      real(8)          :: ztot, zheat
      real(8)          :: axpow(kd)    ! automatic
      real(8)          :: axmid(kd)    ! automatic
      real(8)          :: axmin(kd)    ! automatic
      real(8)          :: axmax(kd)    ! automatic
      integer          :: numax(kd)    ! automatic

      if (kd.le.1) return   ! skip 1D edits

!--- calculate axial averages

      axpow(:)=0.0d0
      axmin(:)=1.0d20
      axmax(:)=0.0d0
      numax(:)=0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
          na=mapcore(ia,ja)
          if (na.eq.0) cycle

          do k=1, kd       ! loop over axial levels
            do j=1, npin
              do i=1, npin
                pp=power(i,j,k,na)
                if (pp.gt.0.0d0) then
                  axpow(k)=axpow(k)+pp
                  numax(k)=numax(k)+1
                  axmin(k)=min(pp,axmin(k))  ! don't count zero power locations
                  axmax(k)=max(pp,axmax(k))
                endif
              enddo
            enddo
          enddo    ! k
        enddo
      enddo

!  Normalize each axial level
!  Note: this will give "odd" results if different axial levels have different number of pins

      pave=0.0d0   ! core average
      zave=0.0d0   ! height of fuel
      nn=0
      do k=1, kd
        if (numax(k).eq.0) then
           axmin(k)=0    ! protect overflow
        else
           axpow(k)=axpow(k)/dble(numax(k))
           pave=pave+axial(k)*axpow(k)
           zave=zave+axial(k)
           if (nn.eq.0) nn=numax(k)       ! save first non-zero level
           if (nn.ne.numax(k)) nn=-1000   ! set error
        endif
      enddo
      if (zave.gt.0.0d0) pave=pave/zave
      if (nn.lt.0) then
        write (*,*) '*** 1D axial distributions not calculated because some ***'
        write (*,*) '***    axial levels have different numbers of pins     ***'
        return
      endif

!--- calculate heated length

      ztot=0.0d0
      zheat=0.0d0
      do k=1, kd
        ztot=ztot+axial(k)
        if (axpow(k).gt.0.0d0) zheat=zheat+axial(k)
      enddo

!--- find axial midpoints for plotting

      z1=0.0d0
      do k=1, kd
        z2=z1+axial(k)
        axmid(k)=0.5d0*(z1+z2)
        z1=z2
      enddo

      write (*,'(/,1x,a,/,1x,a)') trim(title), '1D Axial Distribution'
      write (*,'(a)')      '    K   Elev         Delta      Average     Minimum     Maximum'
      if (pave.lt.0.001d0) then  ! print exponent
        do k=kd, 1, -1
          write (*,210) k, axmid(k), axial(k), axpow(k), axmin(k), axmax(k)
        enddo
      else
        do k=kd, 1, -1
          write (*,200) k, axmid(k), axial(k), axpow(k), axmin(k), axmax(k)
        enddo
      endif
      write (*,205) pave

      write (*,140) 'total ', ztot
      write (*,140) 'heated', zheat
  140 format (2x,a,' axial length ', f12.4,' cm')

  200 format (i5,f10.4,6f12.4)
  210 format (i5,f10.4,f12.4,1p,6e14.5)
  205 format (2x,'ave',22x,f12.4)

      return
      end subroutine print1d

!=======================================================================
!
!  Collapse 3D edits to 2D
!
! 05/19/2017 - took out 3D max/min edits, these are calculated elsewhere
!            - added option to collapse with loadings
!
!=======================================================================
      subroutine collapse2d(label, npin, kd, nassm, axial, power, pow2d, ifload)
      use mod_coregeom, only : pinload
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)
      real(8), intent(out):: pow2d(npin,npin,nassm)
      logical, intent(in) :: ifload   ! flag to weight with loadings

!--- local

      integer :: i, j, k
      integer :: na
      integer :: nhot
      integer :: k2min(3)
      integer :: k2max(3)
      real(8) :: pp, pl
      real(8) :: zave, zrod
      real(8) :: c2min, c2max

      pow2d(:,:,:)=0.0d0

      write (*,'(/,1x,a,a)') 'Collapsing 3D distributions to 2D - ', trim(label)

!--- collapse to 2D  (loading)

      if (ifload) then

        write (*,*) 'performing axial collapse with pin loadings'
        do na=1, nassm     ! loop over assemblies
          do j=1, npin
            do i=1, npin
              zrod=0.0d0
              zave=0.0d0   ! 2D axial height with power
              do k=1, kd   ! loop over axial levels
                pp=power(i,j,k,na)
                if (pp.gt.0.0d0) then
                  pl=pinload(i,j,k,na)
                  zrod=zrod+axial(k)*pl*pp
                  zave=zave+axial(k)*pl
                  if (pl.eq.0.0d0) write (*,*) 'WARNING: zero loading used in collapse'
                endif
              enddo
              if (zave.gt.0.0d0) pow2d(i,j,na)=zrod/zave
            enddo
          enddo
        enddo      ! na

      else

!--- collapse to 2D  (no loading)

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
              if (zave.gt.0.0d0) pow2d(i,j,na)=zrod/zave
            enddo
          enddo
        enddo      ! na

      endif   ! pin loadings

!--- find max/min

      c2max=0.0d0        ! 2D core max
      c2min=1.0d20       ! 2D core min
      k2max(:)=0         ! 2D core max i, j, n
      k2min(:)=0         ! 2D core min i, j, n
      nhot=0

      do na=1, nassm     ! loop over assemblies
        do j=1, npin
          do i=1, npin
            pp=pow2d(i,j,na)
            if (pp.gt.0.0d0) then 
              nhot=nhot+1
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

      if (nhot.eq.0) then    ! protect edits
        c2min=0.0d0
      endif

      write (*,190) 'max', c2max, k2max(:)
      write (*,190) 'min', c2min, k2min(:)
  190 format ('   2D ',a,' in core =', f10.4,' at (i,j,na)  ', 2i4,4x,i4)

      return
      end subroutine collapse2d

!=======================================================================
!
!  Subroutine to print 2D Assembly Maps (assembly average values)
!
!=======================================================================
      subroutine print2d_assm_map(title, npin, nassm, pow2)
      use mod_coregeom, only : icore, jcore, mapcore, xlabel, ylabel
      implicit none
      integer, intent(in) :: npin, nassm
      real(8), intent(in) :: pow2(npin,npin,nassm)
      character(len=*), intent(in) :: title

      character(len=8) :: fmt

      integer          :: i, j, nn
      integer          :: nnsave
      integer          :: ia, ja, na
      real(8)          :: pmin, pmax
      real(8)          :: pp
      real(8)          :: passm(nassm)  ! automatic

      logical  :: ifbw  ! flag if assemblies have different number of pins

      ifbw=.false.

      write (*,'(/,1x,2a,/)') trim(title), ' - 2D assembly average (2RPF)'

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
          ifbw=.true.
!         write (0,*) 'debug: assembly ', na, nn, nnsave
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
      write (*,125) nn
      if (ifbw) then
        write (*,*) 'WARNING: average of assemblies cannot be calculated because assemblies have different number of rods'
      else
        write (*,130) 'average', pp
!!      if (abs(pp-1.0d0).gt.0.0001) write (0,*) '***** check normalization *****'   ! msg only valid for pin powers
      endif
      write (*,130) 'maximum', pmax
      write (*,130) 'minimum', pmin
  125 format (4x,'number of assemblies in full-core ', i0)
  130 format (4x,a,' assembly power  ', f8.4)

      write (*,*)

      fmt='(f8.4)'
      if (pp.gt.100.0d0) fmt='(f8.2)'

!--- print map

      write (*,'("  **  ",50(4x,a2,2x))') (xlabel(i),i=1,icore)   ! labels are 2

      do j=1, jcore
        write (*,'(2x,a2,"- ")',advance='no') ylabel(j)
        do i=1, icore
           if (mapcore(i,j).eq.0) then
             write (*,'(8x)',advance='no')
           else
             write (*,fmt,advance='no') passm(mapcore(i,j))
           endif
        enddo
        write (*,*)
      enddo

      return
      end subroutine print2d_assm_map
!=======================================================================
!
!  Subroutine to print 2D max pin in assembly maps (2PIN)
!
!  2017/05/18 - added for Jim
!
!=======================================================================
      subroutine print2d_2pin_map(title, npin, nassm, power, kd)
      use mod_coregeom, only : icore, jcore, mapcore, axial, xlabel, ylabel
      implicit none
      integer, intent(in) :: npin, nassm
      integer, intent(in) :: kd
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

      character(len=8) :: fmt

      integer          :: i, j, k
      integer          :: ia, ja, na
      real(8)          :: pp
      real(8)          :: zave, zlen
      real(8)          :: pin2d(icore,jcore)  ! automatic
      real(8)          :: core2d      ! core max

      write (*,'(/,1x,2a,/)') trim(title), ' - max 2D rod in each assembly (2PIN)'

!--- calculate max pin in each assembly - use mapcore in case this is qtr-core

      pin2d(:,:)=0.0d0
      core2d=0.0d0

      do ja=1, jcore
        do ia=1, icore
          na=mapcore(ia,ja)
          if (na.gt.0) then

            do j=1, npin
              do i=1, npin
                zave=0.0d0
                zlen=0.0d0
                do k=1, kd   ! loop over axial levels
                  pp=power(i,j,k,na)
                  if (pp.gt.0.0d0) then
                    zave=zave+axial(k)*pp
                    zlen=zlen+axial(k)
                  endif
                enddo
                if (zlen.gt.0.0d0) then
                  zave=zave/zlen
                  pin2d(ia,ja)=max(pin2d(ia,ja),zave)
                endif
              enddo
            enddo

            core2d=max(core2d,pin2d(ia,ja))   ! core max

          endif  ! na
        enddo    ! ia
      enddo      ! ja

      fmt='(f8.4)'
      if (core2d.gt.100.0d0) fmt='(f8.2)'

!--- print map

      write (*,'("  **  ",50(4x,a2,2x))') (xlabel(i),i=1,icore)   ! labels are 2

      do ja=1, jcore
        write (*,'(2x,a2,"- ")',advance='no') ylabel(ja)
        do ia=1, icore
           if (mapcore(ia,ja).eq.0) then
             write (*,'(8x)',advance='no')
           else
             write (*,fmt,advance='no') pin2d(ia,ja)
           endif
        enddo
        write (*,*)
      enddo

      write (*,*) 'WARNING: 2PIN maps do not use loadings if collapsing exposures'

      return
      end subroutine print2d_2pin_map
!=======================================================================
!
!  Subroutine to calculate core-wide axial offset
!
!  General version that considers qtr-symmetry
!
!=======================================================================
      subroutine calc_axoff (title, npin, kd, nassm, power, axoff)
      use mod_coregeom, only : icore, jcore, mapcore, axial
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(out):: axoff  ! Axial offset
      character(len=*), intent(in) :: title

!--- local

      real(8)  :: zmesh(0:kd)    ! automatic

      integer  :: i, j, k
      integer  :: ia, ja, na
      integer  :: kmid
      real(8)  :: chi
      real(8)  :: pp
      real(8)  :: ptop, pbot
      real(8)  :: ztop, zbot
      real(8)  :: zmid
      real(8)  :: zold

      axoff=0.0d0
      if (kd.eq.1) return   ! exit for 2D problems

!--- calculate midpoint

      zold=0.0d0
      do k=1, kd
        zmesh(k)=zold+axial(k)
        zold=zmesh(k)
      enddo
      zmid=zmesh(kd)*0.5d0

!d    write (*,*) 'debug: core height ', zmesh(kd)
!d    write (*,*) 'debug: core mid    ', zmid

      chi=1.0d0    ! set if boundaries align
      kmid=0
      zold=0.0d0
      do k=1, kd
        if (zmid.ge.zold .and. zmid.lt.zmesh(k)) then    ! fix 2019/02/21
          kmid=k
          chi=(zmid-zmesh(k-1))/axial(k)
        endif
        zold=zmesh(k)  ! protect bounds when k=1
      enddo

!d    write (*,*) 'debug: kmid ', kmid
!d    write (*,*) 'debug: chi  ', chi

!--- calculate 3D statistics

      ptop=0.0d0     ! axial offset edits
      pbot=0.0d0
      ztop=0.0d0
      zbot=0.0d0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        do k=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                if (k.lt.kmid) then
                  pbot=pbot+axial(k)*pp
                  zbot=zbot+axial(k)
                elseif (k.gt.kmid) then
                  ptop=ptop+axial(k)*pp
                  ztop=ztop+axial(k)
                else    ! core midplane in node
                  pbot=pbot+axial(k)*chi*pp
                  zbot=zbot+axial(k)*chi
                  ptop=ptop+axial(k)*(1.0d0-chi)*pp
                  ztop=ztop+axial(k)*(1.0d0-chi)
                endif
              endif
            enddo  ! i
          enddo    ! j
        enddo      ! k

      enddo    ! ia
      enddo    ! ja

!d    write (*,*) 'debug: zbot ', zbot
!d    write (*,*) 'debug: ztop ', ztop

      ptop=ptop/ztop
      pbot=pbot/zbot
      axoff=(ptop-pbot)/(ptop+pbot)*100.0d0

      write (*,210) trim(title), axoff
  210 format (2x,a,' A/O  ',f10.4,' %')

      return
      end subroutine calc_axoff
