!=======================================================================
!
!  Subroutines for writing 1D and 2D maps to output
!
!  Copyright (c) 2014-2015 Core Physics, Inc.
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
!=======================================================================
      subroutine stat3d (title, npin, kd, nassm, icore, jcore, mapcore, axial, power, pave, pmax)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      integer, intent(in) :: icore, jcore
      integer, intent(in) :: mapcore(icore,jcore)
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(out):: pave   ! average
      real(8), intent(out):: pmax   ! max
      character(len=*), intent(in) :: title

!--- local

      integer  :: i, j, k
      integer  :: ia, ja, na, kpin
      real(8)  :: pp
      real(8)  :: zlen, zave      ! axial values
      real(8)  :: c3min, c3max    ! 3D values
      integer  :: k3min(4)        ! 3D min locations
      integer  :: k3max(4)        ! 3D max locations

!--- calculate 3D statistics

      pave=0.0d0   ! output average
      pmax=0.0d0   ! output maximum

      kpin=0
      zave=0.0d0
      zlen=0.0d0
      c3min=1.0d20
      c3max=0.0d0
      k3min(:)=0
      k3max(:)=0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        do k=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
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

      pave=zave    ! return max value
      pmax=c3max   ! return max value

!--- check for very small values

      do k=1, kd
        do na=1, nassm
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                if (pp.lt.zave*1.0d-4) then
                   write (*,*) 'small value detected ', i, j, k, na, pp
                endif
              endif
            enddo
          enddo
        enddo
      enddo


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

      do k=1, kd
        do ia=1, nassm
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,ia)
              if (pp.gt.0.0d0) then
                if (pp.lt.zave*1.0d-4) then
                   write (*,*) 'small value detected ', i, j, k, ia, pp
                endif
              endif
            enddo
          enddo
        enddo
      enddo


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
      character(len=150) :: line
      character(len=8)   :: fmt

      if (npin.gt.20) then
        write (*,*) '**** too many pins across to print nice maps ****'
        return
      endif

!--- calculate maximum and set the format correctly

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
      do i=1, npin
        line=' '
        m=0
        do j=1, npin
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
      subroutine print1d(title, npin, kd, nassm, power, axial)
      implicit none
      character(len=*), intent(in) :: title
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)

      integer          :: i, j, k, nn
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

!--- calculate axial averages **** this will be wrong if qtr-core and assemblies repeated

      axpow(:)=0.0d0
      axmin(:)=1.0d20
      axmax(:)=0.0d0
      numax(:)=0

      do na=1, nassm
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
!=======================================================================
      subroutine collapse2d(label, npin, kd, nassm, axial, power, pow2d)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      real(8), intent(in) :: axial(kd)
      real(8), intent(out):: pow2d(npin,npin,nassm)

!--- local

      integer :: i, j, k
      integer :: na
      integer :: k3min(4)
      integer :: k3max(4)
      integer :: k2min(3)
      integer :: k2max(3)
      real(8) :: pp
      real(8) :: zave, zrod
      real(8) :: c3min, c3max
      real(8) :: c2min, c2max
      real(8) :: cave, clen

      pow2d(:,:,:)=0.0d0

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

      write (*,'(/,1x,a,a)') 'Collapsing 3D edits to 2D - ', trim(label)

!--- collapse to 2D

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
            if (zave.gt.0.0d0) then
              pow2d(i,j,na)=zrod/zave
              pp=pow2d(i,j,na)
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

      write (*,180) 'max', c3max, k3max(:)
      write (*,180) 'min', c3min, k3min(:)
      write (*,190) 'max', c2max, k2max(:)
      write (*,190) 'min', c2min, k2min(:)

  180 format (' 3D ',a,' in core =', f10.4,' at (i,j,k,na)', 4i4)
  190 format (' 2D ',a,' in core =', f10.4,' at (i,j,na)  ', 2i4,4x,i4)

      return
      end subroutine collapse2d

