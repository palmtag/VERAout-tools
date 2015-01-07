!=======================================================================
!
!  Subroutines for writing 1D and 2D maps to output
!
!  Copyright (c) 2014 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  Subroutine to print 3D statistics
!
!=======================================================================
      subroutine stat3d (title, npin, kd, nassm, axial, power)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: power(npin,npin,kd,nassm)
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
      end subroutine stat3d

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

      integer            :: i, j, ia
      integer            :: k
      real(8), allocatable :: ptemp(:,:,:)

!--- create 2D map and pass to edit routine

      allocate (ptemp(npin,npin,nassm))
      do ia=1, nassm
        do j=1, npin
          do i=1, npin
            ptemp(i,j,ia)=power(i,j,kd,ia)
          enddo
        enddo
      enddo
     
      write (*,*) 
      write (*,*) 'Exit values level ', kd
      k=1  ! for 2D map
      call print_pin_map(title, npin, k, nassm, ptemp)

      deallocate (ptemp)

      return
      end subroutine print_exit_map
!=======================================================================
!
!  Subroutine to print 3D pin powers by assembly using pretty output
!
!=======================================================================
      subroutine print_pin_map(title, npin, kd, nassm, power)
      implicit none
      integer, intent(in) :: npin, kd, nassm
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer            :: i, j, m
      integer            :: ia, klev, kpin
      real(8)            :: pp, pave, pmin, pmax
      real(8)            :: zmax   ! 3D max
      character(len=150) :: line
      character(len=8)   :: fmt

!--- calculate 3D statistics to set the format correctly
!---  **** average does not use axial weighting

      zmax=0.0d0

      do ia=1, nassm
        do klev=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,klev,ia)
              zmax=max(zmax,pp)
            enddo
          enddo
        enddo
      enddo

!!    write (*,'(/,1x,a)') trim(title)
!!    write (*,'(a,2f12.4)') ' Max ', zmax

!--- write maps


      fmt='(f7.4)'
      if (zmax.gt.  10.0d0) fmt='(f7.3)'
      if (zmax.gt. 100.0d0) fmt='(f7.2)'
      if (zmax.gt.1000.0d0) fmt='(f7.1)'

!d    write (*,*) 'debug: zmax = ', zmax
!d    write (*,*) 'debug: fmt  = ', fmt 

      if (npin.gt.20) then
        write (*,*) '**** too many pins across to print nice maps ****'
        return
      endif

!**** To-do: read labels from input and use in edits


      do ia=1, nassm        ! loop over assemblies
        do klev=kd, 1, -1   ! loop over axial levels
          write (*,'(/,1x,a)') trim(title)
!!        write (*,'(1x,a,i3,2a)') 'Assembly ', ia,'    type ', trim(assmmap(nassm))   ! no types to mpact
          write (*,'(1x,a,i3,2a)') 'Assembly ', ia
          if (kd.gt.1) then
            write (*,'(  1x,a,i3)') 'Level    ', klev
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
!d        write (*,*) 'debug average = ', pave

        enddo    ! klev
      enddo      ! ia
! 210 format (5x,'number of hot pins',i4,'    average=',f7.4,'   min=', f7.4,'   max=', f7.4)
  210 format (5x,'number of hot pins',i4,'    average=',f12.4,'   min=', f12.4,'   max=', f12.4)

      return
      end subroutine print_pin_map
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

